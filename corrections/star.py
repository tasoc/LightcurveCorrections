import numpy as np
"""
A class for storing star metadata for ensemble photometry.

.. codeauthor:: Oliver J. Hall
.. codeauthor:: Fillipe Pereira
"""


class Star(object):
	"""
	A class for storing metadata for ensemble photometry.

	Input:
		time (ndarray): Time data for a stellar flux timeseries.
		flux (ndarray): Flux data for a stellar flux timeseries.
	"""
	def __init__(self, time, flux):
		self.time = time
		self.flux = flux
		self.init_stats()

	def init_stats(self):
		"""
		A function to build star metadata off of the flux timeseries.

		Properties:
			fmean (float): Mean of the flux
			fstd (float): Standard deviation of second differences of the flux
			frange (float): absolute relative 5-95th percentile flux range normalised by mean
			drange (float): relative differenced (whitened) standard deviation
		"""

		# Mean flux
		self.fmean = np.mean(self.flux)

		# Standard deviation of twice-differenced (whitened) time series
		self.fstd = np.std(np.diff(np.diff(self.flux)))

		# Relative 5-95 percentile flux range
		frange = np.percentile(self.flux,95) - np.percentile(self.flux,5)
		frange = frange / np.mean(self.flux)
		self.frange = abs(frange)

		# Relative differenced (whitened) standard deviation
		drange = np.std(np.diff(self.flux)) / np.mean(self.flux)
		self.drange = abs(drange)
