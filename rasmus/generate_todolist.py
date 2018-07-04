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
import sys

if __name__ == '__main__':
	# Make sure some directories exist:
	os.makedirs('data/noisy_by_sectors', exist_ok=True)
	os.makedirs('data/clean_by_sectors', exist_ok=True)

	progress = 0
	if os.path.exists('progress.txt'):
		with open('progress.txt', 'r') as fid:
			progress = int(fid.readline())

	print(progress)

	pri = {}
	starid = 0
	for starlist_file in ('data/Data_Batch_TDA.txt', 'data/Data_Batch_TDA_solar.txt', 'data/Data_Batch_TDA_correction.txt'):

		starlist = np.genfromtxt(starlist_file, delimiter=',', dtype=None)

		for k, star in enumerate(starlist):
			print(star)

			# Don't trust Mikkel!
			#starname = star[0]
			#if not isinstance(starname, six.string_types): starname = starname.decode("utf-8") # For Python 3
			#starid = int(starname[4:])
			starid = starid+1

			if starid <= progress:
				continue

			startype = star[-1]
			if not isinstance(startype, six.string_types): startype = startype.decode("utf-8") # For Python 3

			tmag = star[1]

			if starid <= 15000:
				file_to_load = os.path.join('noisy_files_TDA', 'Star%d.noisy' % starid)
			elif starid <= 17001:
				file_to_load = os.path.join('noisy_files_TDA_solar', 'Star%d.noisy' % (starid-15001))
			elif starid <= 37001:
				file_to_load = os.path.join('noisy_files_TDA_correction', 'Star%d.noisy' % (starid-17001))
			else:
				raise Exception("Here there be dragons!")

			data = np.loadtxt(os.path.join('data', file_to_load))
			data_clean = np.loadtxt(os.path.join('data', file_to_load.replace('noisy', 'clean')))

			sector = np.floor(data[:,0] / 27.4) + 1
			sectors = [int(s) for s in np.unique(sector)]

			# Just because Mikkel can not be trusted:
			#if star[2] == 1800:
			#	datasource = 'ffi'
			#else:
			#	datasource = 'tpf'
			if (data[1,0] - data[0,0])*86400 > 1000:
				datasource = 'ffi'
			else:
				datasource = 'tpf'

			# Extract the camera from the lattitude:
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
				print('Star%d-sector%02d.noisy' % (starid, s))

				indx = (sector == s)
				data_sector = data[indx, :]

				# Save files cut up into sectors:
				lightcurve = 'Star%d-sector%02d.noisy' % (starid, s)
				np.savetxt(os.path.join('data', 'noisy_by_sectors', lightcurve), data_sector)
				np.savetxt(os.path.join('data', 'clean_by_sectors', 'Star%d-sector%02d.clean' % (starid, s)), data_clean[indx, :])

				sqlite_file = 'todo-sector%02d.sqlite' % s
				if not os.path.exists(sqlite_file):
					conn = sqlite3.connect(sqlite_file)
					cursor = conn.cursor()

					cursor.execute("""CREATE TABLE todolist (
						priority BIGINT PRIMARY KEY NOT NULL,
						starid BIGINT NOT NULL,
						datasource TEXT NOT NULL DEFAULT 'ffi',
						camera INT NOT NULL,
						ccd INT NOT NULL,
						method TEXT DEFAULT NULL,
						tmag REAL,
						status INT DEFAULT NULL,
						cbv_area INT NOT NULL,
						lightcurve TEXT DEFAULT NULL
					);""")

					cursor.execute("""CREATE TABLE diagnostics (
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
					cursor.execute("CREATE UNIQUE INDEX priority_idx ON todolist (priority);")
					cursor.execute("CREATE INDEX starid_datasource_idx ON todolist (starid, datasource);") # FIXME: Should be "UNIQUE", but something is weird in ETE-6?!
					cursor.execute("CREATE INDEX status_idx ON todolist (status);")
					cursor.execute("CREATE INDEX starid_idx ON todolist (starid);")
					cursor.execute("CREATE INDEX variability_idx ON diagnostics (variability);")
					conn.commit()

					pri[s] = 0
				else:
					conn = sqlite3.connect(sqlite_file)
					cursor = conn.cursor()
					if s not in pri:
						cursor.execute("SELECT MAX(priority) FROM todolist;")
						pri[s] = int(cursor.fetchone()[0])

				pri[s] += 1
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

				cursor.execute("INSERT INTO todolist (priority,starid,tmag,datasource,status,camera,ccd,cbv_area,lightcurve) VALUES (?,?,?,?,0,?,?,?,?);", (
					pri[s],
					starid,
					tmag,
					datasource,
					camera,
					1,
					cbv_area,
					lightcurve
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
				cursor.close()
				conn.close()

			# Write progress to file:
			with open('progress.txt', 'w') as fid:
				fid.write(str(starid))

	print("DONE.")
