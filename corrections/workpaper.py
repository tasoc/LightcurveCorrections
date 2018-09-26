"""
WORKPAPER - TESS Corrections

quick and easy end-to-end implementation of a single batch of lightcurves
through the general detrending strategy combining filters and CBV output.
"""

import os


#------------------------------------------------------------------------------
def fileIntake():
    # generalized method for taking in lightcurves from various sources
    # step 1: identify the sources
    # step 2: filter on source type 
    # step 3: read in relevant data
    pass

#------------------------------------------------------------------------------
def fitsIntake():
    # more specialized method for taking in lightcurves from FITS sources
    #
    #
    raise NotImplementedError

#------------------------------------------------------------------------------
def saveOutput():
    # generalized method for saving various corrections related output files 
    # step 1: identify type of output(s)
    # step 2: filter on source type and call out to helpers as needed
    # step 3: save and verify
    pass

#------------------------------------------------------------------------------
def parseCommandline():
    # generalized method for defining and parsing input commands; better Base*?
    # NOTE: most "run" or "prepare" codes use parsers, so create a general list
    #       ('debug', 'quiet', etc.) like this and add to it case-by-case
    # step 1: add arguments to parser
    # step 2: check argument validity 
    # step 3: setup based on arguments - logging, etc.
    pass

#------------------------------------------------------------------------------
