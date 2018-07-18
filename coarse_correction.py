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
        #
        ###

        return STATUS.OK
