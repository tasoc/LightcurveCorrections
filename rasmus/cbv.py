#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

from __future__ import division, with_statement, print_function, absolute_import
from six.moves import range
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import os
import glob
from sklearn.decomposition import PCA
from bottleneck import allnan, nansum, move_median, nanmedian, nanstd
from scipy.optimize import minimize
from scipy.stats import pearsonr
from scipy.interpolate import pchip_interpolate
import itertools
#from photometry.utilities import move_median_central
from tqdm import tqdm

#Hello
#Hello again

#------------------------------------------------------------------------------
def _move_median_central_1d(x, width_points):
	y = move_median(x, width_points, min_count=1)
	y = np.roll(y, -width_points//2+1)
	for k in range(width_points//2+1):
		y[k] = nanmedian(x[:(k+2)])
		y[-(k+1)] = nanmedian(x[-(k+2):])
	return y

#------------------------------------------------------------------------------
def move_median_central(x, width_points, axis=0):
	return np.apply_along_axis(_move_median_central_1d, axis, x, width_points)

#------------------------------------------------------------------------------
def pearson(x, y):
	indx = np.isfinite(x) & np.isfinite(y)
	r, _ = pearsonr(x[indx], y[indx]) # Second output (p-value) is not used
	return r

#------------------------------------------------------------------------------
class CBV(object):

	def __init__(self, filepath):
		self.cbv = np.load(filepath)

	def mdl(self, coeffs):
		coeffs = np.atleast_1d(coeffs)
		m = np.empty(self.cbv.shape[0], dtype='float64')
		m.fill(coeffs[-1])
		for k in range(len(coeffs)-1):
			m += coeffs[k] * self.cbv[:, k]
		return m

	def _lhood(self, coeffs, flux):
		return nansum((flux - self.mdl(coeffs))**2)

	def fit(self, flux, Ncbvs=2, sigma_clip=4.0, maxiter=3):

		# Find the median flux, as it is used for
		# initial guesses later on:
		median_flux = nanmedian(flux)

		# Start looping over the number of CBVs to include:
		bic = np.empty(self.cbv.shape[1]+1, dtype='float64')
		solutions = []
		for Ncbvs in range(self.cbv.shape[1]+1):

			# Initial guesses for coefficients:
			coeffs0 = np.zeros(Ncbvs + 1, dtype='float64')
			coeffs0[-1] = median_flux

			iters = 0
			fluxi = np.copy(flux)
			while iters <= maxiter:
				iters += 1

				# Do the fit:
				res = minimize(self._lhood, coeffs0, args=(fluxi, ), method='Powell')
				flux_filter = self.mdl(res.x)

				# Do robust sigma clipping:
				absdev = np.abs( fluxi - flux_filter )
				mad = 1.4826*nanmedian(absdev)
				indx = np.greater(absdev, sigma_clip*mad, where=np.isfinite(fluxi))
				if np.any(indx):
					fluxi[indx] = np.nan
				else:
					break

			# Calculate the Bayesian Information Criterion (BIC) and store the solution:
			bic[Ncbvs] = np.log(np.sum(np.isfinite(fluxi)))*len(coeffs0) + res.fun
			solutions.append(res)

		# Use the solution which minimizes the BIC:
		indx = np.argmin(bic)
		flux_filter = self.mdl(solutions[indx].x)

		#plt.figure()
		#plt.plot(bic, '.-')
		#plt.show()

		return flux_filter

#------------------------------------------------------------------------------
if __name__ == '__main__':

	# Pick a sector, any sector....
	sector = 2

	# Other settings:
	threshold_variability = 1.3
	threshold_correlation = 0.5

	# Remove old plots:
	os.makedirs("plots/sector%02d/" % sector, exist_ok=True)
	for f in glob.iglob("plots/sector%02d/*.png" % sector):
		os.remove(f)

	# Open the TODO file for that sector:
	filepath_todo = 'todo-sector%02d.sqlite' % sector
	conn = sqlite3.connect(filepath_todo)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()

	# Get list of CBV areas:
	cursor.execute("SELECT DISTINCT cbv_area FROM todolist ORDER BY cbv_area;")
	cbv_areas = [int(row[0]) for row in cursor.fetchall()]
	print(cbv_areas)

	# Loop through the CBV areas:
	# - or run them in parallel - whatever you like!
	for cbv_area in cbv_areas:

		#---------------------------------------------------------------------------------------------------------
		# CALCULATE CBV FOR THIS CBV-AREA
		#---------------------------------------------------------------------------------------------------------

		print("We are running CBV_AREA=%d" % cbv_area)
		camera = np.floor(cbv_area/100)
		#print(camera)

		tmpfile = 'mat-sector%02d-%d.npz' % (sector, cbv_area)
		if os.path.exists(tmpfile):
			print("Loading existing file...")
			data = np.load(tmpfile)
			mat = data['mat']
			priorities = data['priorities']
			stds = data['stds']

		else:
			# Find the median of the variabilities:
			# SQLite does not have a median function so we are going to
			# load all the values into an array and make Python do the
			# heavy lifting.
			cursor.execute("""SELECT variability FROM todolist LEFT JOIN diagnostics ON todolist.priority=diagnostics.priority WHERE
				datasource='ffi'
				AND status=0
				AND CAST(cbv_area/100 AS INT)=?;""", (camera, ))
			variability = np.array([row[0] for row in cursor.fetchall()], dtype='float64')
			median_variability = nanmedian(variability)

			# Plot the distribution of variability for all stars:
			fig = plt.figure()
			ax = fig.add_subplot(111)
			ax.hist(variability/median_variability, bins=np.logspace(np.log10(0.1), np.log10(1000.0), 50))
			ax.axvline(threshold_variability, color='r')
			ax.set_xscale('log')
			ax.set_xlabel('Variability')
			fig.savefig('plots/sector%02d/variability-area%d.png' % (sector, cbv_area))
			plt.close(fig)

			# Get the list of star that we are going to load in the lightcurves for:
			# We have put a hard limit of a maximum of 20000 targets for now.
			cursor.execute("""SELECT todolist.starid,todolist.priority,mean_flux,variance FROM todolist LEFT JOIN diagnostics ON todolist.priority=diagnostics.priority WHERE
				datasource='ffi'
				AND status=0
				AND CAST(cbv_area/100 AS INT)=?
				AND variability < ?
			ORDER BY variability ASC LIMIT 30000;""", (camera, threshold_variability*median_variability))
			stars = cursor.fetchall()

			# Number of stars returned:
			Nstars = len(stars)

			# Load the very first timeseries only to find the number of timestamps.
			# When using FITS files in the future, this could simply be loaded from the header.
			time = np.loadtxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (stars[0]['starid'], sector)), usecols=(0, ), unpack=True)
			Ntimes = len(time)

			print("Matrix size: %d x %d" % (Nstars, Ntimes))

			# Make the matrix that will hold all the lightcurves:
			print("Loading in lightcurves...")
			mat = np.empty((Nstars, Ntimes), dtype='float64')
			mat.fill(np.nan)
			stds = np.empty(Nstars, dtype='float64')
			priorities = np.empty(Nstars, dtype='int64')
			for k, star in tqdm(enumerate(stars), total=Nstars):
				priorities[k] = star['priority']
				starid = star['starid']

				flux = np.loadtxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, sector)), usecols=(1, ), unpack=True)
				#print(flux.shape)

				# Normalize the data and store it in the rows of the matrix:
				mat[k, :] = flux / star['mean_flux'] - 1.0
				stds[k] = np.sqrt(star['variance'])

			# Only start calculating correlations if we are actually filtering using them:
			if threshold_correlation < 1.0:
				file_correlations = 'correlations-sector%02d-%d.npy' % (sector, cbv_area)
				if os.path.exists(file_correlations):
					correlations = np.load(file_correlations)
				else:
					# Calculate the correlation matrix between all lightcurves:
					print("Calculating correlations...")
					correlations = np.empty((Nstars, Nstars), dtype='float64')
					np.fill_diagonal(correlations, np.nan) # Put NaNs on the diagonal
					for i, j in tqdm(itertools.combinations(range(Nstars), 2), total=0.5*Nstars**2-Nstars):
						r = pearson(mat[i, :]/stds[i], mat[j, :]/stds[j])
						correlations[i,j] = correlations[j,i] = np.abs(r)

					np.save(file_correlations, correlations)

				# Find the median absolute correlation between each lightcurve and all other lightcurves:
				c = nanmedian(correlations, axis=0)

				# Indicies that would sort the lightcurves by correlations in descending order:
				indx = np.argsort(c)[::-1]
				indx = indx[:int(threshold_correlation*Nstars)]

				# Only keep the top 50% of the lightcurves that are most correlated:
				priorities = priorities[indx]
				mat = mat[indx, :]

				# Clean up a bit:
				del correlations, c, indx

			# Save something for debugging:
			np.savez('mat-sector%02d-%d.npz' % (sector, cbv_area), mat=mat, priorities=priorities, stds=stds)

		# Print the final shape of the matrix:
		print("Matrix size: %d x %d" % mat.shape)

		# Simple low-pass filter of the individual targets:
		#mat = move_median_central(mat, 48, axis=1)

		# Find columns where all stars have NaNs and remove them:
		indx_nancol = allnan(mat, axis=0)
		Ntimes = mat.shape[1]
		mat = mat[:, ~indx_nancol]

		cadenceno = np.arange(mat.shape[1])

		# TODO: Is this even needed? Or should it be done earlier?
		print("Gap-filling lightcurves...")
		for k in tqdm(range(mat.shape[0]), total=mat.shape[0]):

			mat[k, :] /= stds[k]

			# Fill out missing values by interpolating the lightcurve:
			indx = np.isfinite(mat[k, :])
			mat[k, ~indx] = pchip_interpolate(cadenceno[indx], mat[k, indx], cadenceno[~indx])

		#replace(mat, np.nan, 0) # Replace NaNs with zero... We should do something better...

		# Calculate the principle components:
		print("Doing Principle Component Analysis...")
		pca = PCA(n_components=8)
		pca.fit(mat)

		# Not very clever, but puts NaNs back into the CBVs:
		# For some reason I also choose to transpose the CBV matrix
		cbv = np.empty((Ntimes, 8), dtype='float64')
		cbv.fill(np.nan)
		cbv[~indx_nancol, :] = np.transpose(pca.components_)

		# Plot the "effectiveness" of each CBV:
		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.plot(np.arange(1, cbv.shape[1]+1), pca.explained_variance_ratio_, '.-')
		ax.set_xlabel('CBV number')
		ax.set_ylabel('Variance explained ratio')
		fig.savefig('plots/sector%02d/cbv-expvar-area%d.png' % (sector, cbv_area))
		plt.close(fig)

		# Plot all the CBVs:
		fig, axes = plt.subplots(4, 2, figsize=(12, 8))
		for k, ax in enumerate(axes.flatten()):
			ax.plot(cbv[:, k], '-')
			ax.set_title('Basis Vector %d' % (k+1))
		plt.tight_layout()
		fig.savefig('plots/sector%02d/cbvs-area%d.png' % (sector, cbv_area))
		plt.close(fig)

		# Save the CBV to file:
		np.save('cbv-sector%02d-%d.npy' % (sector, cbv_area), cbv)

		#---------------------------------------------------------------------------------------------------------
		# CORRECTING STARS
		#---------------------------------------------------------------------------------------------------------

		print("CORRECTING STARS...")

		# Query for all stars, no matter what variability and so on
		cursor.execute("""SELECT starid FROM todolist WHERE
			datasource='ffi'
			AND cbv_area=?
			AND status=0;""", (cbv_area, ))
		stars = cursor.fetchall()

		# Load the cbv from file:
		cbv = CBV('cbv-sector%02d-%d.npy' % (sector, cbv_area))

		for star in stars:
			starid = star['starid']


			time, flux = np.loadtxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, sector)), usecols=(0, 1), unpack=True)
			time_clean, flux_clean = np.loadtxt(os.path.join('data', 'clean_by_sectors', 'Star%d-sector%02d.clean' % (starid, sector)), usecols=(0, 1), unpack=True)

			# Fit the CBV to the flux:
			flux_filter = cbv.fit(flux, Ncbvs=4)

			fig = plt.figure()
			ax1 = fig.add_subplot(211)
			ax1.plot(time, flux)
			ax1.plot(time, flux_filter)
			ax1.set_xticks([])
			ax2 = fig.add_subplot(212)
			ax2.plot(time, flux/flux_filter-1)
			ax2.plot(time_clean, flux_clean/nanmedian(flux_clean)-1, alpha=0.5)
			ax2.set_xlabel('Time')
			plt.tight_layout()
			fig.savefig('plots/sector%02d/star%d.png' % (sector, starid))
			plt.close(fig)
