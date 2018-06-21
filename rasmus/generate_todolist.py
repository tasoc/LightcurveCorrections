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
	starlist = np.genfromtxt('data/Data_Batch_TDA.txt', delimiter=',', dtype=None)

	pri = {}
	for k, star in enumerate(starlist):
		print(star)

		#starname = star[0]
		#if not isinstance(starname, six.string_types): starname = starname.decode("utf-8") # For Python 3
		#starid = int(starname[4:])
		starid = k+1

		tmag = star[1]
		if star[2] == 1800:
			datasource = 'ffi'
		else:
			datasource = 'tpf'

		data = np.loadtxt(os.path.join('data', 'noisy_files_TDA', 'Star%d.noisy' % starid))
		sector = np.floor(data[:,0] / 27.4) + 1
		sectors = [int(s) for s in np.unique(sector)]

		for s in sectors:
			if s not in pri: pri[s] = 1
			indx = (sector == s)
			data_sector = data[indx, :]
			print('Star%d-sector%02d.noisy' % (starid, s))

			np.savetxt(os.path.join('data', 'noisy_by_sectors', 'Star%d-sector%02d.noisy' % (starid, s)), data_sector)

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
				mask_size INT,
				pos_row REAL,
				pos_column REAL,
				contamination REAL,
				stamp_resizes INT,
				errors TEXT
			);""")
			conn.commit()

			mean_flux = nanmedian(data_sector[:,1])
			variance = nansum((data_sector[:,1] - mean_flux)**2) / (data_sector.shape[0] - 1)

			elaptime = np.random.normal(3.14, 0.5)
			pos_row = np.random.uniform(-0.5, 2047.5)
			pos_column = np.random.uniform(-0.5, 2047.5)
			Npixels = np.interp(tmag, np.array([8.0, 9.0, 10.0, 12.0, 14.0, 16.0]), np.array([350.0, 200.0, 125.0, 100.0, 50.0, 40.0]))

			if pos_row > 1024 and pos_column > 1024:
				cbv_area = 114
			elif pos_row > 1024 and pos_column <= 1024:
				cbv_area = 113
			elif pos_row <= 1024 and pos_column > 1024:
				cbv_area = 112
			else:
				cbv_area = 111

			#if mean_flux < 0:
			#	plt.figure()
			#	plt.plot(data[:,0], data[:,1])
			#	plt.show()

			cursor.execute("INSERT INTO todolist (priority,starid,tmag,datasource,status,camera,ccd,cbv_area) VALUES ({priority:d},{starid:d},{tmag:f},'{datasource:s}',0,{camera:d},{ccd:d},{cbv_area:d});".format(
				priority=pri[s],
				starid=starid,
				tmag=tmag,
				datasource=datasource,
				camera=1,
				ccd=1,
				cbv_area=cbv_area
			))
			cursor.execute("INSERT INTO diagnostics (priority,starid,elaptime,mean_flux,variance,mask_size,pos_row,pos_column,contamination,stamp_resizes) VALUES ({priority:d},{starid:d},{elaptime:f},{mean_flux:f},{variance:f},{mask_size:d},{pos_row:f},{pos_column:f},0.0,0);".format(
				priority=pri[s],
				starid=starid,
				elaptime=elaptime,
				mean_flux=mean_flux,
				variance=variance,
				mask_size=int(Npixels),
				pos_row=pos_row,
				pos_column=pos_column
			))

			conn.commit()

			pri[s] += 1

	cursor.close()
	conn.close()