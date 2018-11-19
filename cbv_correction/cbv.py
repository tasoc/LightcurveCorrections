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
import sys
import glob
from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score
from bottleneck import allnan, nansum, move_median, nanmedian, nanstd
from scipy.optimize import minimize
from scipy.stats import pearsonr, entropy
from scipy.interpolate import pchip_interpolate
import itertools
from statsmodels.nonparametric.kde import KDEUnivariate as KDE
from scipy.special import xlogy

import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module="scipy.stats") # they are simply annoying!

#from photometry.utilities import move_median_central
from tqdm import tqdm
plt.ioff()
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
def compute_scores(X, n_components):
    pca = PCA(svd_solver='full')

    pca_scores = []
    for n in n_components:
        pca.n_components = n
        pca_scores.append(np.mean(cross_val_score(pca, X, cv=5)))

    return pca_scores

#------------------------------------------------------------------------------
def rms(x, **kwargs):
	return np.sqrt(nansum(x**2, **kwargs)/len(x))

#------------------------------------------------------------------------------
def hpd(data, level) :
    """ The Highest Posterior Density (credible) interval of data at level level.
  
    :param data: sequence of real values
    :param level: (0 < level < 1)
    """ 
    
    d = list(data)
    d.sort()
  
    nData = len(data)
    nIn = int(round(level * nData))
    if nIn < 2 :
      raise RuntimeError("not enough data")
    
    i = 0
    r = d[i+nIn-1] - d[i]
    for k in range(len(d) - (nIn - 1)) :
      rk = d[k+nIn-1] - d[k]
      if rk < r :
        r = rk
        i = k
  
    assert 0 <= i <= i+nIn-1 < len(d)
    
    return (d[i], d[i+nIn-1], i, i+nIn-1)
    

#------------------------------------------------------------------------------
def compute_entopy(U):

	HGauss0 = 0.5 + 0.5*np.log(2*np.pi)
		
	nSingVals = U.shape[1]
	H = np.empty(nSingVals, dtype='float64')
	
#	fig, axes = plt.subplots(4, 2, figsize=(12, 8))
	
	
#	for iBasisVector, ax in enumerate(axes.flatten()):
	for iBasisVector in range(nSingVals):
		
#		print(np.sum(np.abs(U[:,iBasisVector]>0)), U.shape[0])
#		pdf, xedges = np.histogram(np.abs(U[:, iBasisVector]), bins=nBins, density=True)
#		x = xedges[1:] - xedges[:-1]
		
		kde = KDE(np.abs(U[:, iBasisVector]))
		kde.fit(gridsize=1000)
		
		pdf = kde.density
		x = kde.support
		
		dx = x[1]-x[0]
		# Calculate the Gaussian entropy
		pdfMean = nansum(x * pdf)*dx
#		pdfMean = kde.support[np.argmax(kde.density)]
		sigma = np.sqrt( nansum(((x-pdfMean)**2) * pdf) * dx )
		
#		HPD = hpd(np.abs(U[:, iBasisVector]), 0.682689492137)
#		print(HPD)
#		sigma = (HPD[1] - HPD[0])/2
		HGauss = HGauss0 + np.log(sigma)

		# Calculate vMatrix entropy
		pdf_pos = (pdf>0)
#		HVMatrix = entropy(pdf[pdf_pos] ) * dx 
		HVMatrix = -np.sum(xlogy(pdf[pdf_pos], pdf[pdf_pos])) * dx

		# Returned entropy is difference between V-Matrix entropy and Gaussian entropy of similar width (sigma)
		H[iBasisVector] = HVMatrix - HGauss	
		
#		print(pdfMean, sigma, HGauss, HVMatrix)
#		ax.plot(kde.support, kde.density)
#	plt.show()	
	return H
	

#------------------------------------------------------------------------------
def clean_cbv(Matrix, n_components, ent_limit=-1.5, targ_limit=50):
	
	
	# Calculate the principle components:
	print("Doing Principle Component Analysis...")
	pca = PCA(n_components)
	U, _, _ = pca._fit(Matrix)
	
	
	
	Ent = compute_entopy(U)
	print('Entropy start:', Ent)
	
	targets_removed = 0
	components = np.arange(n_components)
		
	while np.any(Ent<ent_limit):
		com = components[(Ent<ent_limit)][0]
		
		# Remove highest relative weight target
		m = nanmedian(U[:, com])
		s = 1.46*nanmedian(np.abs(U[:, com] - m))
		dev = np.abs(U[:, com] - m) / s

		idx0 = np.argmax(dev)
		
		# Remove high weight target
#		idx0 = np.argmax(np.abs(U[:,com]))	
		star_no = np.ones(U.shape[0], dtype=bool)
		star_no[idx0] = False
		print('removing star ', idx0)
		
		Matrix = Matrix[star_no, :]
		U, _, _ = pca._fit(Matrix)

		targets_removed += 1
		
		if targets_removed>targ_limit:
			break
		
		Ent = compute_entopy(U)
		print('Entropy:', Ent)
		
	
	
	print('Targets removed ', targets_removed)
	return Matrix




def lc_matrix(sector, cbv_area):
	
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
			AND status=1
			AND cbv_area=?;""", (cbv_area, ))
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
			AND status=1
			AND cbv_area=?
			AND variability < ?
		ORDER BY variability ASC LIMIT 30000;""", (cbv_area, threshold_variability*median_variability))
		stars = cursor.fetchall()

		# Number of stars returned:
		Nstars = len(stars)

		# Load the very first timeseries only to find the number of timestamps.
		# When using FITS files in the future, this could simply be loaded from the header.
		time = np.loadtxt(os.path.join('sysnoise', 'Star%d.sysnoise' % (stars[0]['starid'],)), usecols=(0, ), unpack=True)
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

			flux = np.loadtxt(os.path.join('sysnoise', 'Star%d.sysnoise' % (starid,)), usecols=(1, ), unpack=True)
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
			#TODO: remove based on threshold value? rather than just % of stars

			# Only keep the top 50% of the lightcurves that are most correlated:
			priorities = priorities[indx]
			mat = mat[indx, :]

			# Clean up a bit:
			del correlations, c, indx

		# Save something for debugging:
		np.savez('mat-sector%02d-%d.npz' % (sector, cbv_area), mat=mat, priorities=priorities, stds=stds)

	return mat, priorities, stds
	
def lc_matrix_clean(sector, cbv_area):
	
	
	print('Running matrix clean')
	tmpfile = 'mat-sector%02d-%d_clean.npz' % (sector, cbv_area)
	if os.path.exists(tmpfile):
		print("Loading existing file...")
		data = np.load(tmpfile)
		mat = data['mat']
		priorities = data['priorities']
		stds = data['stds']
		
		Ntimes = data['Ntimes']
		indx_nancol = data['indx_nancol']

	else:
		# Compute light curve correlation matrix
		mat0, priorities, stds = lc_matrix(sector, cbv_area)
		
		# Print the final shape of the matrix:
		print("Matrix size: %d x %d" % mat0.shape)

		# Simple low-pass filter of the individual targets:
		#mat = move_median_central(mat, 48, axis=1)

		# Find columns where all stars have NaNs and remove them:
		indx_nancol = allnan(mat0, axis=0)
		Ntimes = mat0.shape[1]
		mat = mat0[:, ~indx_nancol]

		cadenceno = np.arange(mat.shape[1])

		# TODO: Is this even needed? Or should it be done earlier?
		print("Gap-filling lightcurves...")
		for k in tqdm(range(mat.shape[0]), total=mat.shape[0]):

			mat[k, :] /= stds[k]

			# Fill out missing values by interpolating the lightcurve:
			indx = np.isfinite(mat[k, :])
			# Do inpainting??
			mat[k, ~indx] = pchip_interpolate(cadenceno[indx], mat[k, indx], cadenceno[~indx])

		# Save something for debugging:
		np.savez('mat-sector%02d-%d_clean.npz' % (sector, cbv_area), mat=mat, priorities=priorities, stds=stds, indx_nancol=indx_nancol, Ntimes=Ntimes)

	return mat, priorities, stds, indx_nancol, Ntimes
	
#------------------------------------------------------------------------------
class CBV(object):

	def __init__(self, filepath):
		self.cbv = np.load(filepath)
		
	def lsfit(self, flux):
		
		idx = np.isfinite(self.cbv[:,0]) & np.isfinite(flux)
		A = self.cbv[idx,:]
		coeffs = np.linalg.lstsq(A, flux[idx], rcond=None)[0]
		
		return coeffs
		
	def mdl(self, coeffs):
		coeffs = np.atleast_1d(coeffs)
#		m = np.empty(self.cbv.shape[0], dtype='float64')
#		m.fill(coeffs[-1])
		m = np.ones(self.cbv.shape[0], dtype='float64')
		
#		for k in range(len(coeffs)-1):
		for k in range(len(coeffs)):
			m += coeffs[k] * self.cbv[:, k]
		return m
	

	def _lhood(self, coeffs, flux):
		return nansum((flux - self.mdl(coeffs))**2)
	

	def fit(self, flux, Numcbvs=2, sigma_clip=4.0, maxiter=3, use_bic=True):

		# Find the median flux to normalise light curve
		median_flux = nanmedian(flux)
		
		if Numcbvs is None:
			Numcbvs = self.cbv.shape[1]   
			
		if use_bic:	
			# Start looping over the number of CBVs to include:
			bic = np.empty(Numcbvs+1, dtype='float64')
			solutions = []
			
			# Test a range of CBVs from 0 to Numcbvs
			Nstart = 0 
		else:
			# Test only fit with Numcbvs
			Nstart = Numcbvs
	
	
	
		for Ncbvs in range(Nstart, Numcbvs+1):

			# Initial guesses for coefficients:
#			coeffs0 = np.zeros(Ncbvs + 1, dtype='float64')
			coeffs0 = np.zeros(Ncbvs, dtype='float64')
			coeffs0[0] = 0.3
#			coeffs0[-1] = 1

			iters = 0
			fluxi = np.copy(flux) / median_flux 
			while iters <= maxiter:
				iters += 1

				# Do the fit:
				res = minimize(self._lhood, coeffs0, args=(fluxi, ), method='Powell')
				flux_filter = self.mdl(res.x)

#				res = self.lsfit(fluxi)
#				flux_filter = self.mdl(res)

				# Do robust sigma clipping:
				absdev = np.abs( fluxi - flux_filter )
				mad = 1.4826*nanmedian(absdev)
				indx = np.greater(absdev, sigma_clip*mad, where=np.isfinite(fluxi))
				
				# Update guess for next iteration
#				coeffs0 = res.x
				
				if np.any(indx):
					fluxi[indx] = np.nan
				else:
					break

			if use_bic:
				# Calculate the Bayesian Information Criterion (BIC) and store the solution:
				bic[Ncbvs] = np.log(np.sum(np.isfinite(fluxi)))*len(coeffs0) + res.fun
				solutions.append(res)
			
							

		if use_bic:
			# Use the solution which minimizes the BIC:
			indx = np.argmin(bic)
			res_final = solutions[indx].x
			flux_filter = self.mdl(res_final)  * median_flux

		else:
			res_final = res.x
#			res_final = res
			flux_filter = self.mdl(res_final)  * median_flux

		#plt.figure()
		#plt.plot(bic, '.-')
		#plt.show()

		return flux_filter, res_final

#------------------------------------------------------------------------------
if __name__ == '__main__':

	import pandas as pd
	# Pick a sector, any sector....
	sector = 0
	n_components = 8

	# Other settings:
	threshold_variability = 1.3
	threshold_correlation = 0.5
	threshold_snrtest = 5.0

	
	# Open the TODO file for that sector:
	filepath_todo = 'todo.sqlite'
	conn = sqlite3.connect(filepath_todo)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()

	# Get list of CBV areas:
	cursor.execute("SELECT DISTINCT cbv_area FROM todolist ORDER BY cbv_area;")
	cbv_areas = [int(row[0]) for row in cursor.fetchall()]
	print(cbv_areas)

	# Loop through the CBV areas:
	# - or run them in parallel - whatever you like!
	for ii, cbv_area in enumerate(cbv_areas):
		
		n_components = 8
		
		if ii<4:
			continue
		
		# Remove old plots:
		if os.path.exists("plots/sector%02d/" % sector):
			for f in glob.iglob("plots/sector%02d/*%d.png" % (sector, cbv_area)):
				os.remove(f)		
		else:    
			os.makedirs("plots/sector%02d/" % sector)
			
		# Remove old plots:
		if os.path.exists("plots/sector%02d/stars_area%d/" % (sector, cbv_area)):
			for f in glob.iglob("plots/sector%02d/stars_area%d/*.png" % (sector, cbv_area)):
				os.remove(f)		
		else:    
			os.makedirs("plots/sector%02d/stars_area%d/" % (sector, cbv_area))



		# Extract or compute cleaned and gapgilled light curve matrix
		mat0, priorities, stds, indx_nancol, Ntimes = lc_matrix_clean(sector, cbv_area)
		
		# Calculate initial CBVs
		pca0 = PCA(n_components)
		U0, _, _ = pca0._fit(mat0)
		# Not very clever, but puts NaNs back into the CBVs:
		# For some reason I also choose to transpose the CBV matrix
		cbv0 = np.empty((Ntimes, n_components), dtype='float64')
		cbv0.fill(np.nan)
		cbv0[~indx_nancol, :] = np.transpose(pca0.components_)
		
		
		print('Cleaning matrix for CBV - remove single dominant contributions')
		# Clean away targets that contribute significantly as a single star to a given CBV (based on entropy)
		mat = clean_cbv(mat0, n_components, ent_limit=-2, targ_limit=150)
		

		# Calculate the principle components of cleaned matrix
		print("Doing Principle Component Analysis...")
		pca = PCA(n_components)
		U, _, _ = pca._fit(mat)
		
		
		# Not very clever, but puts NaNs back into the CBVs:
		# For some reason I also choose to transpose the CBV matrix
		cbv1 = np.empty((Ntimes, n_components), dtype='float64')
		cbv1.fill(np.nan)
		cbv1[~indx_nancol, :] = np.transpose(pca.components_)
		
		
		# Signal-to-Noise test:
		A_signal = rms(cbv1, axis=0)
		A_noise = rms(np.diff(cbv1, axis=0), axis=0)
		snr = 10 * np.log10( A_signal**2 / A_noise**2 )
		indx_lowsnr = (snr < threshold_snrtest)
		if np.any(indx_lowsnr):
			print("Rejecting %d CBVs based on SNR test" % np.sum(indx_lowsnr))
			cbv = cbv1[:, ~indx_lowsnr]
		else:
			cbv = cbv1
			
		# Update maximum number of components	
		n_components = cbv.shape[1]
		print('New max number of components: ', n_components)
		
		# Save the CBV to file:
		np.save('cbv-sector%02d-%d.npy' % (sector, cbv_area), cbv)
		
		max_components=20
		n_cbv_components = np.arange(max_components, dtype=int)
		pca_scores = compute_scores(mat, n_cbv_components)
		
		# Plot the "effectiveness" of each CBV:
		fig0 = plt.figure(figsize=(12,8))
		ax0 = fig0.add_subplot(121)
		ax02 = fig0.add_subplot(122)
		ax0.plot(n_cbv_components, pca_scores, 'b', label='PCA scores')
		ax0.set_xlabel('nb of components')
		ax0.set_ylabel('CV scores')
		ax0.legend(loc='lower right')
		
		ax02.plot(np.arange(1, cbv0.shape[1]+1), pca.explained_variance_ratio_, '.-')
		ax02.axvline(x=cbv.shape[1]+0.5, ls='--', color='k')
		ax02.set_xlabel('CBV number')
		ax02.set_ylabel('Variance explained ratio')
		
		fig0.savefig('plots/sector%02d/cbv-perf-area%d.png' % (sector, cbv_area))
		plt.close(fig0)
		

		# Plot all the CBVs:
		fig, axes = plt.subplots(4, 2, figsize=(12, 8))
		fig2, axes2 = plt.subplots(4, 2, figsize=(12, 8))
		for k, ax in enumerate(axes.flatten()):
			ax.plot(cbv0[:, k], 'r-')		
			if indx_lowsnr[k]:
				col = 'c'
			else:
				col = 'k'
			ax.plot(cbv1[:, k], ls='-', color=col)	
			ax.set_title('Basis Vector %d' % (k+1))
			
			
		for k, ax in enumerate(axes2.flatten()):	
			ax.plot(-np.abs(U0[:, k]), 'r-')
			ax.plot(np.abs(U[:, k]), 'k-')
			ax.set_title('Basis Vector %d' % (k+1))
		plt.tight_layout()
		fig.savefig('plots/sector%02d/cbvs-area%d.png' % (sector, cbv_area))
		fig2.savefig('plots/sector%02d/U_cbvs-area%d.png' % (sector, cbv_area))
		plt.close(fig)

		


#		sys.exit()
		#---------------------------------------------------------------------------------------------------------
		# CORRECTING STARS
		#---------------------------------------------------------------------------------------------------------

		print("CORRECTING STARS...")

		# Query for all stars, no matter what variability and so on
		cursor.execute("""SELECT starid FROM todolist WHERE
			datasource='ffi'
			AND cbv_area=?
			AND status=1;""", (cbv_area, ))
		stars = cursor.fetchall()

		# Load the cbv from file:
		cbv = CBV('cbv-sector%02d-%d.npy' % (sector, cbv_area))
		
		results = np.zeros([len(stars), n_components+2])

		for kk, star in enumerate(stars):
			starid = star['starid']
			print(starid, kk, len(stars))

			data = pd.read_csv(os.path.join('sysnoise', 'Star%d.sysnoise' % (starid,)),  usecols=(0, 1), skiprows=6, sep=' ', header=None, names=['Time', 'Flux'])
			time, flux = data['Time'].values, data['Flux'].values

			# Fit the CBV to the flux:
			flux_filter, res = cbv.fit(flux, Numcbvs=n_components, use_bic=False)
			
			res = np.array([res,]).flatten()
			results[kk, 0] = starid
			results[kk, 1:len(res)+1] = res
			print(res)
			
			
			# Plot comparison between clean and corrected data
			data_clean = pd.read_csv(os.path.join('noisy', 'Star%d.noisy' % (starid,)),  usecols=(0, 1), skiprows=6, sep=' ', header=None, names=['Time', 'Flux'])
			time_clean, flux_clean = data_clean['Time'].values, data_clean['Flux'].values

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
			fig.savefig('plots/sector%02d/stars_area%d/star%d.png' % (sector, cbv_area, starid))
			plt.close('all')




		np.savez('mat-sector%02d-%d_free_weights.npz' % (sector, cbv_area), res=results)


		fig = plt.figure(figsize=(15,6))
		ax = fig.add_subplot(121)
		ax2 = fig.add_subplot(122)
		for kk in range(1,n_components+1):
			idx = np.nonzero(results[:, kk])
			r = results[idx, kk]
			idx2 = (r>np.percentile(r, 10)) & (r<np.percentile(r, 90))
			kde = KDE(r[idx2])
			kde.fit(gridsize=5000)
			
			ax.plot(kde.support*1e5, kde.density/np.max(kde.density), label='CBV ' + str(kk))
			
			err = nanmedian(np.abs(r[idx2] - nanmedian(r[idx2]))) * 1e5
			ax2.errorbar(kk, kde.support[np.argmax(kde.density)]*1e5, yerr=err, marker='o', color='k')
		ax.set_xlabel('CBV weight')
		ax2.set_ylabel('CBV weight')
		ax2.set_xlabel('CBV')
		ax.legend()
		fig.savefig('plots/sector%02d/weights-sector%02d-%d.png' % (sector, sector, cbv_area))
		plt.close('all')

#		sys.exit()


