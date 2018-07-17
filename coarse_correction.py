#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
TODO: important info goes here

Structure based on `photometry.py` by Rasmus Handberg <rasmush@phys.au.dk>

"""

import numpy as np 
import logging

class CoarseCorrection(BaseCorrection):

    def __init__(self, *args, **kwargs):
        # NOTE: placeholder
        super().__init__()

        # create instances? one for cotrend, one for detrend, etc.?
        self.cotrend = CLC(self.camera, self.ccd)

    def do_correction(self):
        """ Coarse Correction """ 

        logger = logging.getLogger(__name__)

        # Generate list of stars to fit:
        cat = self.catalog
        # TODO: look into this. I think it makes sense to make a list of everyone nearby
        #       on area in order to do other ensemble work in detrending (instead of cotrending CBVs)

        # TODO: make sure the target is in the catalog itself! if it's being used that way?
        #       otherwise it should be in self.etc but might be relevant later??

        ###
        # NOTE: Placeholder for coarse classification logic
        #
        # There may be an argument for using a filter (e.g., on highly variable targets)
        # here and running corrections differently even on the coarse level.
        # That filter probably belongs here, where we can direct the order/choice of 
        # correction calls below.
        #
        # Example: 
        # if method is None and self.variability < threshold:
        #    do_cotrend()
        #    do_detrend()
        # elif self.variability > threshold:
        #    do_detrend()   # no cotrend because the CBV's don't improve EB's, w/e
        #
        ###

        # NOTE: Placeholder
        # TODO: Call the funtions!

        return STATUS.OK

    def do_cotrend():
        """
        Use the CBVs to fit against the given target and remove systematics 

        TODO: implement this
        """

        #
        raise NotImplementedError
    
    def do_detrend():
        """
        Use Derek's detrending fit 

        TODO: implement this
        TODO: decide how many stars to load in for the ensemble attempt - do try something along the lines of:
            "SELECT * FROM todolist LEFT JOIN diagnostics ON todolist.priority=diagnositics.priority WHERE ((list of parameters goes here))"

            At least one of those parameters needs to be a sphere_distance() from the target that's being detrended, and
            we probably want a list of at least 100 in case any don't match with the later checks
        """

        #