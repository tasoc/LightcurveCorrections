# -*- coding: utf-8 -*-
"""
Command-line utility to run TESS detrend correction from command-line.

Note:
    This file is intended as an addition to the code
    `tess_detrend` by Derek Buzasi <dbuzasi@fgcu.edu>

    It allows for small tests of single target inputs from
    a variety of formats, including FITS and text files,
    and provides an option for debugging and maintenance
Structure inspired by `tessphot` by Rasmus Handberg <rasmush@phys.au.dk>

"""

# photometry runs under Python 2, so compatibility is an issue without this
from __future__ import with_statement, print_function
import os
import argparse
import fnmatch
import logging
import functools

#------------------------------------------------------------------------------
if __name__ == '__main__':
    # list of files to run through
    fitsfiles = [] 

    # allowing for the option of running on a single star/folder 
    # command line parameters (input, output, format, ...)
    parser = argparse.ArgumentParser(description='Run detrending correction for a single target from the TESS Photometry pipeline')
    parser.add_argument('-i', '--input', help='input folder, default local: ' + os.getcwd(), default=os.getcwd())
    parser.add_argument('-o', '--output', help='output path, default local: ' + os.getcwd(), default=os.getcwd())
    parser.add_argument('-d', '--debug', help='Print debug messages.', action='store_true')
    parser.add_argument('-q', '--quiet', help='Only report warnings and errors.', action='store_true')
    parser.add_argument('-f', '--format', help='use alternate output format, default tab-delimited', nargs='?', const='latex', default='tab')	
            # NOTE: probably don't want a latex default, but it's an okay placeholder for now
    args = parser.parse_args()

    # Make sure appropriate settings are supplied, or that defaults are acceptable (print/log)

    # Set logging level - keeps consistency with logging in run_tessphot
    logging_level = logging.INFO 
    if args.quiet:
        logging_level = logging.WARNING 
    elif args.debug:
        logging_level = logging.DEBUG 

    # Setup logging 
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger = logging.getLogger(__name__)
    logger.addHandler(console)
    logger.setLevel(logging_level)
    logger_parent = logging.getLogger('correction')
    logger_parent.addHandler(console)
    logger_parent.setLevel(logging_level)

    # run_tessphot uses environment variables to set in/out folders, can use for defaults
    test_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), 'input'))
    if args.input:
        input_foler = test_folder
    else:
        input_folder = os.environ.get('TESSCORR_INPUT', test_folder)
    output_folder = os.environ.get('TESSCORR_OUTPUT', os.path.abspath('.'))
    logger.info("Loading input data from '%s'", input_folder)
    logger.info("Putting output data in '%s'", output_folder)

    # "Create partial function of tessphot, setting the common keywords:"
	# f = functools.partial(tessphot, input_folder=input_folder, output_folder=output_folder, plot=args.plot)
    ###### TODO: figure out what that means - how does functools.partial() work? ########

    # Run the program 
    with TaskManager(input_folder) as tm:
        # NOTE: photometry can run on single targets, some detrending can too; photometry's option is:
        #  if args.starid is not None:
        #      task = tm.get_task(starid=args.starid)
        # would need some elif/else statement to follow
        task = tm.get_task()
        del task['priority']
        # corr = f(**task) # NOTE: see TODO above

    # Write out results?
	if not args.quiet:
		print("=======================")
		print("STATUS: {0}".format(pho.status.name))
