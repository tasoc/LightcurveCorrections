
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import os
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from bottleneck import replace, allnan, nansum, move_median, nanmedian
from scipy.optimize import minimize
from scipy.stats import pearsonr
import itertools
#from photometry.utilities import move_median_central

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
	r, pval = pearsonr(x[indx], y[indx])
	return r

#------------------------------------------------------------------------------
class CBV(object):

	def __init__(self, filepath):
		self.cbv = np.load(filepath)

	def mdl(self, coeffs):
		m = coeffs[-1]
		for k in range(len(coeffs)-1):
			m += coeffs[k]*self.cbv[:,k]
		return m

	def _lhood(self, coeffs, flux):
		return nansum((flux - self.mdl(coeffs))**2)

	def fit(self, flux, Ncbvs=2):
		
		coeffs0 = np.zeros(Ncbvs + 1)
		res = minimize(self._lhood, coeffs0, args=(flux, ), method='Powell')

		flux_filter = self.mdl(res.x)

		return flux_filter

#------------------------------------------------------------------------------
if __name__ == '__main__':

	# Pick a sector, any sector....
	sector = 1

	# Open the TODO file for that sector:
	filepath_todo = 'todo-sector%02d.sqlite' % sector
	conn = sqlite3.connect(filepath_todo)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()
	
	# Get list of CBV areas:
	cursor.execute("SELECT DISTINCT cbv_area FROM todolist ORDER BY cbv_area;")
	cbv_areas = [int(row[0]) for row in cursor.fetchall()]
	print(cbv_areas)
	colors = ['r', 'b', 'g', 'y']

	# Loop through the CBV areas:
	# - or run them in parallel - whatever you like!
	for cbv_area in cbv_areas:

		#---------------------------------------------------------------------------------------------------------
		# CALCULATE CBV FOR THIS CBV-AREA
		#---------------------------------------------------------------------------------------------------------

		print("We are running CBV_AREA=%d" % cbv_area)

		cursor.execute("""SELECT * FROM todolist LEFT JOIN diagnostics ON todolist.priority=diagnostics.priority WHERE
			datasource='ffi'
			AND cbv_area=?
			AND status=0;""", (cbv_area, ))
		stars = cursor.fetchall()
		
		Nstars = len(stars)
		Ntimes = 1316 # We should maybe store this somewhere???
		print(Nstars)

		# Make the matrix that will hold all the lightcurves:
		mat = np.empty((Nstars, Ntimes), dtype='float64')
		mat.fill(np.nan)

		variability = np.zeros(Nstars)
		for k, star in enumerate(stars):
			starid = star['starid']

			data = np.loadtxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, sector)))
			#print(data.shape)

			# This could be done in the photometry code as well:
			flux = data[:,1] / star['mean_flux']
			indx = np.isfinite(flux)
			p = np.polyfit(data[indx,0], flux[indx], 3)
			variability[k] = np.nanstd(flux - np.polyval(p, data[:,0]))

			# Normalize the data and store it in the rows of the matrix:
			mat[k, :] = data[:,1] / star['mean_flux']

		# 
		variability = variability/np.median(variability)

		fig = plt.figure()
		ax = fig.add_subplot(111)
		ax.hist(variability, bins=np.logspace(np.log10(0.1), np.log10(1000.0), 50))
		ax.set_xscale('log')
		ax.set_xlabel('Variability')
		plt.show()

		# Filter out stars that are variable:
		indx_quiet = (variability < 1.5)
		mat = mat[indx_quiet, :]

		# Calculate the correlation matrix between all lightcurves:
		N = mat.shape[0]
		correlations = np.empty((N, N), dtype='float64')
		np.fill_diagonal(correlations, np.nan) # Put NaNs on the diagonal
		for i, j in itertools.combinations(range(N), 2):
			r = pearson(mat[i, :], mat[j, :])
			correlations[i,j] = correlations[j,i] = np.abs(r)

		# Find the median absolute correlation between each lightcurve and all other lightcurves:
		c = nanmedian(correlations, axis=0)
		
		# Indicies that would sort the lightcurves by correlations in descending order:
		indx = np.argsort(c)[::-1]

		# Only keep the top 50% of the lightcurves that are most correlated:
		mat = mat[indx[:int(0.75*N)], :]

		# Print the final shape of the matrix:
		print(mat.shape)

		# Simple low-pass filter of the individual targets:
		mat = move_median_central(mat, 48, axis=1)	

		# Replace NaNs with zero... We should do something better...
		indx_nancol = allnan(mat, axis=0) 
		mat = mat[:, ~indx_nancol]
		replace(mat, np.nan, 0)
		
		# Calculate the principle components:
		#mat = StandardScaler(copy=False).fit_transform(mat)
		pca = PCA(n_components=5)
		pca.fit(mat)
		cbv = pca.components_
		cbv = np.transpose(cbv)

		# Not very clever, but puts NaNs back into the CBVs:
		mat = np.empty((Ntimes, 5), dtype='float64')
		mat.fill(np.nan)
		mat[~indx_nancol, :] = cbv
		cbv = mat

		# Save the CBV to file:
		np.save('cbv-%d.npy' % cbv_area, cbv)

		#---------------------------------------------------------------------------------------------------------
		# CORRECTING STARS
		#---------------------------------------------------------------------------------------------------------

		cursor.execute("""SELECT * FROM todolist WHERE
			datasource='ffi'
			AND cbv_area=?
			AND status=0;""", (cbv_area, ))
		stars = cursor.fetchall()

		# Load the cbv from file:
		cbv = CBV('cbv-%d.npy' % cbv_area)

		for star in stars:
			starid = star['starid']
			time, flux = np.loadtxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, sector)), usecols=(0, 1), unpack=True)
			
			# Fit the CBV to the flux:
			flux_filter = cbv.fit(flux, Ncbvs=5)
			
			fig = plt.figure()
			plt.plot(time, flux)
			plt.plot(time, flux_filter)
			fig.savefig('plots/star%d-sector%02d.png' % (starid, sector))
			plt.close(fig)
	
	plt.show()


