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

        logger.info('STARID = %d', self.starid) #TODO: others

        self._status = STATUS.UNKNOWN
        self._details = {}
        #TODO: others

        # Directory(s) to save output files:
        #TODO: logic for directing output? CBV vs LC
        self.output_folder = os.path.join(
            output_folder,
            # self.datasource, #TODO: if that is included/relevant 
            '{0:011d}'.format(self.starid)[0:5]
        )

        #TODO: Set directory for diagnostics? photometry `plot` option

        # Init table that will be filled with lightcurve stuff:
        self.lightcurve = Table()

        #TODO: This is where photometry differentiates between input datasources (ffi vs tpf); should prob be where we filter on fits vs txt
        # e.g., if self.datasource == 'fits': ------; if self.datasource == 'txt': ------;
        # consider using 'pho' or similar instead of 'fits' --> indicates output from photometry

        #TODO: whitespace/scope verification after filter implementation
        if camera is None or ccd is None:
            #TODO: should probably exist inside of a datasource filter
            raise ValueError("CAMERA and CCD keywords must be provided for %s targets", self.datasource)

        self.camera = camera # TESS camera
        self.ccd = ccd # TESS CCD
        # TODO: redirect to proper path - might be "area" and not "camera_ccd" based? Talk to Rasmus!!
        self.catalog_file = os.path.join(input_folder, 'catalog_camera{0:d}_ccd{1:d}.sqlite'.format(self.camera, self.ccd))

        logger.debug('CAMERA = %s', self.camera)
        logger.debug('CCD = %s', self.ccd)
        logger.debug('Catalog file: %s', self.catalog_file)

        # Load info on target
        conn = sqlite3.connect(self.catalog_file)
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM catalog WHERE starid={0:d};".format(self.starid))
        target = cursor.fetchone()

        if target is None:
            raise IOError("Target could not be found in catalog: {0:d}".format(self.starid))
        self.target_tmag = target['tmag'] # TESS Magnitude of target
        self.target_pos_ra = target['ra']
        self.target_pos_dec = target['decl'] # TODO: check DB - should be just 'dec'?
        # NOTE: neither of these are J2000 - do we need? is it in the corrections todo?

        # TODO: other database info loading
        cursor.close()
        conn.close()


    @property
    def status(self):
        """The status of the corrections. From :py:class:`STATUS`."""
        return self._status

    @property
    def catalog(self):
        """
        Catalog of targets on current CCD area

        The table contains the following columns:
        * ``starid``:
        * ``tmag``:
        * ``ra``:
        * ``dec``:
        * ``ccdarea``:
        * ``variability``:
        * TODO: other info for corrections?

        Returns:
            ``astropy.table.Table``: Table with all the targets on the current CCD area

        Example:
            if ``corr`` is an instance of :py:class:`BaseCorrection`:
                >>> corr.catalog['tmag']
                >>> corr.catalog[('starid', 'tmag', 'row', 'column')]
        """
        
        if not self._catalog:
            dbconn = sqlite3.connect(self.catalog_file)

    def do_corrections(self):
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

        # run the corrections
        self._status = self.do_corrections(*args, **kwargs)

        #check that the status has been changed
        if self._status == STATUS.UNKNOWN:
            raise Exception("STATUS was not set by do_corrections")

        if self._status in (STATUS.OK, STATUS.WARNING):
            self._details[] = self.lightcurve[]
            self._details[] = self.lightcurve[]
            self._details[] = self.lightcurve[]
            # other outputs, including some conditions