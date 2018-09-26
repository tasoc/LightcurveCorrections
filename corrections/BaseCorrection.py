#!/usr/bin/env python
# -*- coding utf-8 -*-
"""
The basic correction class for the TASOC Photomety pipeline
All other specific correction classes will inherit from BaseCorrection.
Structure from `BasePhotometry` by Rasmus Handberg <rasmush@phys.au.dk>

"""

#----- Imports go here --------------------------------------------------------
import logging
import os
import sqlite3

@enum.unique
class STATUS(enum.Enum):
    """ 
    Status indicator of the status of the correction.
    Differentiates between LC and CBV for the calling function.
    """
    UNKNOWN = 0
    STARTED = 9
    OK = 1
    ERROR = 2
    WARNING = 3
    ABORT = 4
    SKIPPED = 5
    #TODO: need some way of telling the caller that the status means a CBV was created 

class BaseCorrection(object):
    """
    The basic correction class for the TASOC Photometry pipeline.
    All other specific correction classes will inherit from BaseCorrection

    Attributes:
        #-------Placeholder-----------------
    """

    def __init__(self, starid, input_folder, output_folder, 
        camera=None, ccd=None):
        #TODO: add other required arguements
        """
        Initialize the correcion object
        
        Parameters:
            starid (int): TIC number of star to be processed
            input_folder (string): Root directory where files are loaded from.
            output_folder (string): Root directory where output files are saved.
            #TODO:------Placeholder------- (non-optional here)
            #FIXME: are camera, ccd really optional? Is there any case where that is not true? Maybe possible to get from `starid`?
            camera (integer, optional): TESS camera (1-4) target observed on (used for CBV area identification)
            ccd (integer, optional): TESS CCD (1-4) target observed on (used for CBV area identification)
            #TODO:------Placeholder------- (optional here)

        Raises:
            IOError: If starid could not be found (TODO: probably other places as well)
            ValueError: (TODO: on a lot of places)
            NotImplementedError: Everywhere a function has a TODO/FIXME tag preventing execution
        """
        
        logger = logging.getLogger(__name__)
        
        # Store the input:
        self.starid = starid
        self.input_folder = input_folder
        #TODO: others


    @property
    def status(self):
        """The status of the corrections. From :py:class:`STATUS`."""
        return self._status

    def do_correction(self):
        # Following photometry structure 
        """
        Apply corrections to Lightcurve.

        This should fill the following
        * self.cotrended
        * self.detrended

        Returns:
            The status of the corrections.
        
        Raises:
            NotImplementedError
        """
        raise NotImplementedError("TODO: A helpful error message goes here")

    def correction():
        """
        Run correction. 

        Will run the :py:func:`do_corrections` method and check 
        some of the output, calculating various performance metrics

        See also:
            :py:func:`do_corrections`
        """

        # run the correction
        self._status = self.do_correction(*args, **kwargs)

        #check that the status has been changed
        if self._status == STATUS.UNKNOWN:
            raise Exception("STATUS was not set by do_corrections")

        if self._status in (STATUS.OK, STATUS.WARNING):
            self._details[] = self.lightcurve[]
            self._details[] = self.lightcurve[]
            self._details[] = self.lightcurve[]
            # other outputs, including some conditions