#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Rasmus Handberg <rasmush@phys.au.dk>
"""

from __future__ import division, with_statement, print_function, absolute_import
import sqlite3
import six
import numpy as np
import os
from bottleneck import nanmedian, nansum
import matplotlib.pyplot as plt

if __name__ == '__main__':
	pri = {}
	for starlist_file in ('data/Data_Batch_TDA.txt', 'data/Data_Batch_TDA_solar.txt'):

		starlist = np.genfromtxt(starlist_file, delimiter=',', dtype=None)

		for k, star in enumerate(starlist):
			print(star)

			# # Don't trust Mikkel!
			#starname = star[0] 
			#if not isinstance(starname, six.string_types): starname = starname.decode("utf-8") # For Python 3
			#starid = int(starname[4:])
			starid = k+1

			startype = star[-1]
			if not isinstance(startype, six.string_types): startype = startype.decode("utf-8") # For Python 3

			tmag = star[1]
			# Don't trust Mikkel!
			#if star[2] == 1800:
			#	datasource = 'ffi'
			#else:
			#	datasource = 'tpf'

			if starid <= 15000:
				file_to_load = os.path.join('noisy_files_TDA', 'Star%d.noisy' % starid)
			else:
				file_to_load = os.path.join('noisy_files_TDA_solar', 'Star%d.noisy' % (starid-15001))
			print(file_to_load)

			data = np.loadtxt(os.path.join('data', file_to_load))
			data_clean = np.loadtxt(os.path.join('data', file_to_load.replace('noisy', 'clean')))

			sector = np.floor(data[:,0] / 27.4) + 1
			sectors = [int(s) for s in np.unique(sector)]

			# Just because Mikkel can not be trusted:
			print((data[1,0] - data[0,0])*86400)
			if (data[1,0] - data[0,0])*86400 > 1000:
				datasource = 'ffi'
			else:
				datasource = 'tpf'

			lat = star[4]
			if lat < 6+24:
				camera = 1
			elif lat < 6+2*24:
				camera = 2
			elif lat < 6+3*24:
				camera = 3
			else:
				camera = 4

			for s in sectors:
				if s not in pri: pri[s] = 1
				indx = (sector == s)
				data_sector = data[indx, :]
				print('Star%d-sector%02d.noisy' % (starid, s))

				np.savetxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, s)), data_sector)
				np.savetxt(os.path.join('data', 'clean_by_sectors', 'Star%d-sector%02d.clean' % (starid, s)), data_clean[indx, :])

				conn = sqlite3.connect('todo-sector%02d.sqlite' % s)
				cursor = conn.cursor()

				cursor.execute("""CREATE TABLE IF NOT EXISTS todolist (
					priority BIGINT PRIMARY KEY NOT NULL,
					starid BIGINT NOT NULL,
					datasource TEXT NOT NULL DEFAULT 'ffi',
					camera INT NOT NULL,
					ccd INT NOT NULL,
					method TEXT DEFAULT NULL,
					tmag REAL,
					status INT DEFAULT NULL,
					cbv_area INT NOT NULL
				);""")

				cursor.execute("""CREATE TABLE IF NOT EXISTS diagnostics (
					priority BIGINT PRIMARY KEY NOT NULL,
					starid BIGINT NOT NULL,
					elaptime REAL NOT NULL,
					mean_flux DOUBLE PRECISION,
					variance DOUBLE PRECISION,
					variability DOUBLE PRECISION,
					mask_size INT,
					pos_row REAL,
					pos_column REAL,
					contamination REAL,
					stamp_resizes INT,
					errors TEXT
				);""")

				# Create the same indicies as is available in the real todolists:
				cursor.execute("CREATE UNIQUE INDEX IF NOT EXISTS priority_idx ON todolist (priority);")
				cursor.execute("CREATE INDEX IF NOT EXISTS starid_datasource_idx ON todolist (starid, datasource);") # FIXME: Should be "UNIQUE", but something is weird in ETE-6?!
				cursor.execute("CREATE INDEX IF NOT EXISTS status_idx ON todolist (status);")
				cursor.execute("CREATE INDEX IF NOT EXISTS starid_idx ON todolist (starid);")
				cursor.execute("CREATE INDEX IF NOT EXISTS variability_idx ON diagnostics (variability);")
				conn.commit()

				mean_flux = nanmedian(data_sector[:,1])
				variance = nansum((data_sector[:,1] - mean_flux)**2) / (data_sector.shape[0] - 1)

				# This could be done in the photometry code as well:
				time = data_sector[:,0]
				flux = data_sector[:,1] / mean_flux
				indx = np.isfinite(flux)
				p = np.polyfit(time[indx], flux[indx], 3)
				variability = np.nanstd(flux - np.polyval(p, time))

				elaptime = np.random.normal(3.14, 0.5)
				pos_row = np.random.uniform(-0.5, 2047.5)
				pos_column = np.random.uniform(-0.5, 2047.5)
				Npixels = np.interp(tmag, np.array([8.0, 9.0, 10.0, 12.0, 14.0, 16.0]), np.array([350.0, 200.0, 125.0, 100.0, 50.0, 40.0]))

				cbv_area = 100*camera
				if pos_row > 1024 and pos_column > 1024:
					cbv_area += 14
				elif pos_row > 1024 and pos_column <= 1024:
					cbv_area += 13
				elif pos_row <= 1024 and pos_column > 1024:
					cbv_area += 12
				else:
					cbv_area += 11

				#if mean_flux < 0:
				#	plt.figure()
				#	plt.plot(data[:,0], data[:,1])
				#	plt.show()

				cursor.execute("INSERT INTO todolist (priority,starid,tmag,datasource,status,camera,ccd,cbv_area) VALUES (?,?,?,?,0,?,?,?);", (
					pri[s],
					starid,
					tmag,
					datasource,
					camera,
					1,
					cbv_area
				))
				cursor.execute("INSERT INTO diagnostics (priority,starid,elaptime,mean_flux,variance,variability,mask_size,pos_row,pos_column,contamination,stamp_resizes) VALUES (?,?,?,?,?,?,?,?,?,0.0,0);", (
					pri[s],
					starid,
					elaptime,
					mean_flux,
					variance,
					variability,
					int(Npixels),
					pos_row,
					pos_column
				))

				conn.commit()

				pri[s] += 1

				cursor.close()
				conn.close()