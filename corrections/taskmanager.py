#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
A TaskManager which keeps track of which targets to process.

.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
.. codeauthor:: Lindsey Carboneau <lmcarboneau@gmail.com>
"""

from __future__ import division, with_statement, print_function, absolute_import
import os
import sqlite3
import logging
import numpy as np 

class TaskManager(object):
	"""
	A TaskManager which keeps track of which targets to process.
	"""

	def __init__(self, todo_file, cleanup=False):
        """
		Initialize the TaskManager which keeps track of which targets to process.

		Parameters:
			todo_file (string): Path to the TODO-file.
			cleanup (boolean): Perform cleanup/optimization of TODO-file before
			                   during initialization. Default=False.

		Raises:
			IOError: If TODO-file could not be found.
		"""

        if os.path.isdir(todo_file):
			todo_file = os.path.join(todo_file, 'todo.sqlite')

		if not os.path.exists(todo_file):
			raise IOError('Could not find TODO-file')

        # NOTE: Placeholder for detrending options list

		# Setup logging:
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
		console = logging.StreamHandler()
		console.setFormatter(formatter)
		self.logger = logging.getLogger(__name__)
		self.logger.addHandler(console)
		self.logger.setLevel(logging.INFO)

		# Load the SQLite file:
		self.conn = sqlite3.connect(todo_file)
		self.conn.row_factory = sqlite3.Row
		self.cursor = self.conn.cursor()

        # Reset the status of everything for a new run:
		# TODO: This should obviously be removed once we start running for real
		self.cursor.execute("UPDATE todolist SET status=NULL;")
		self.cursor.execute("DROP TABLE IF EXISTS diagnostics;")
		self.conn.commit()

		# Create table for diagnostics:
        # TODO: This is a placeholder diagnostics table from the photometry module - update for detrender!!
		self.cursor.execute("""CREATE TABLE IF NOT EXISTS diagnostics (
			priority INT PRIMARY KEY NOT NULL,
			starid BIGINT NOT NULL,
			lightcurve TEXT,
			elaptime REAL NOT NULL,
			mean_flux DOUBLE PRECISION,
			variance DOUBLE PRECISION,
			mask_size INT,
			pos_row REAL,
			pos_column REAL,
			contamination REAL,
			stamp_resizes INT,
			errors TEXT
		);""")
		self.conn.commit()

        # Run a cleanup/optimization of the database before we get started:
		if cleanup:
			self.logger.info("Cleaning TODOLIST before run...")
			try:
				self.conn.isolation_level = None
				self.cursor.execute("VACUUM;")
			except:
				raise
			finally:
				self.conn.isolation_level = ''

	def close(self):
		"""Close TaskManager and all associated objects."""
		self.cursor.close()
		self.conn.close()
		#self.write_summary() # NOTE: photometry writes a summary, right now classifier doesn't

	def __exit__(self, *args):
		self.close()

	def __enter__(self):
		return self
