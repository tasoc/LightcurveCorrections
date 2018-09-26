# -*- coding: utf-8 -*-
"""
TESS Correction - tesscorr.py
Here's the place to put the code to make Derek's code
run nicely alongside the photometry code

following convention:
tessphot.py == TESS photometry --> tesscorr.py == TESS Correction (detrending)
Structure inspired by `tessphot` by Rasmus Handberg <rasmush@phys.au.dk>

"""

from __future__ import with_statement, print_function, absolute_import
import logging
import traceback
# from . import STATUS, etc...

#------------------------------------------------------------------------------
class _CorrErrorDummy(object):
	def __init__(self, *args, **kwargs):
		self.status = STATUS.ERROR
		self._details = {}

#------------------------------------------------------------------------------
def _try_correction(CorrClass, *args, **kwargs):
	logger = logging.getLogger(__name__)
	# try/except for doing correction
	try:
		with CorrClass(*args, **kwargs) as corr:
			corr.correction()

			if corr.status in (STATUS.OK, STATUS.WARNING):
				corr.save_lightcurve()

	except (KeyboardInterrupt, SystemExit):
		logger.info("Stopped by user or system")
		try:
			corr._status = STATUS.ABORT
		except:
			pass
	except:
		logger.exception("Something happened")
		tb = traceback.format_exc().strip()
		try:
			corr._status = STATUS.ERROR
			corr.report_details(error=tb)
		except:
			pass

    try:
		return corr
	except UnboundLocalError:
		return _CorrErrorDummy(*args, **kwargs)
#------------------------------------------------------------------------------
def tesscorr(method=None, *args, **kwargs):
	"""
	Run the detrending correction for a single star and/or create CBVs

	This function will run the specified correction on a given star, 
	creating CBVs for a given area or CCD if they do not already exist.

	Parameters:
	    #NOTE: placeholder
	
	Returns: 
		#NOTE: placeholder
	"""

	logger = logging.getLogger(__name__)

	if method is None:
		# assume general, coarse correction
		corr = _try_correction(CoarseCorrection, *args, **kwargs)

	elif method == 'classification':
		# call out to a 'ClassifiedCorrection' class that will handle the various special cases?
		pass

	elif method == 'final':
		# call out to a 'FinalCorrection' class - may be part of the 'ClassifiedCorrection'?
		# "final" correction may be defined by multiple passes, etc.; may include classification and special handling
		pass

	else:
		raise ValueError("Invalid method: '{0}'".format(method))