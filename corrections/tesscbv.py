# -*- coding: utf-8 -*-
"""
TESS CBV Creation - tesscbv.py
Here's the place to put the code to make 
Cotrending Basis Vectors

following convention:
tessphot.py == TESS photometry --> tesscbv.py == TESS CBV (Cotrending Basis Vectors)
Structure inspired by `tessphot` by Rasmus Handberg <rasmush@phys.au.dk>

.. code author:: Lindsey Carboneau
"""

from __future__ import with_statement, print_function, absolute_import
import logging
import traceback
# from . import STATUS, etc...

#------------------------------------------------------------------------------
class _CBVErrorDummy(object):
	def __init__(self, *args, **kwargs):
		self.status = STATUS.ERROR
		self._details = {}

#------------------------------------------------------------------------------
def _try_cbv(*args, **kwargs):
	logger = logging.getLogger(__name__)
	# try/except for doing correction
	try:
        cbv()
    except:
        #TODO: don't do this
        pass
#------------------------------------------------------------------------------
def tesscbv(method=None, *args, **kwargs):
	"""
	Run the CBV creator

	This function will create CBVs for a given area or CCD if they do not already exist.

	Parameters:
	    #NOTE: placeholder
	
	Returns: 
		#NOTE: placeholder
	"""

	logger = logging.getLogger(__name__)

	if method is None:
		# assume CCD
        pass
	elif method == 'area':
		# do it for an area
		pass

	else:
		raise ValueError("Invalid method: '{0}'".format(method))