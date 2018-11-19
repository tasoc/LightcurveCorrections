#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
.. codeauthor:: Mikkel N. Lund <mikkelnl@phys.au.dk>
"""

from __future__ import division, with_statement, print_function, absolute_import
from six.moves import range
import numpy as np
import sqlite3
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import sys
import glob
from sklearn.decomposition import PCA
from sklearn.model_selection import cross_val_score
from bottleneck import allnan, nansum, move_median, nanmedian, nanstd
from scipy.optimize import minimize
from scipy.stats import pearsonr, entropy
from scipy.interpolate import pchip_interpolate
import itertools
from statsmodels.nonparametric.kde import KDEUnivariate as KDE
from scipy.special import xlogy
import pandas as pd

import warnings
warnings.filterwarnings('ignore', category=FutureWarning, module="scipy.stats") # they are simply annoying!

from tqdm import tqdm

from cbv_util import compute_entopy, _move_median_central_1d, move_median_central, compute_scores, rms
plt.ioff()

# =============================================================================
# 
# =============================================================================

#------------------------------------------------------------------------------
if __name__ == '__main__':
	# Pick a sector, any sector....
	sector = 0
	n_components = 8
	
	# Open the TODO file for that sector:
	filepath_todo = 'todo.sqlite'
	conn = sqlite3.connect(filepath_todo)
	conn.row_factory = sqlite3.Row
	cursor = conn.cursor()
	
	# Get list of CBV areas:
	cursor.execute("SELECT DISTINCT cbv_area FROM todolist ORDER BY cbv_area;")
	cbv_areas = [int(row[0]) for row in cursor.fetchall()]
	print(cbv_areas)
	
	colormap = plt.cm.PuOr #or any other colormap
	normalize1 = colors.Normalize(vmin=0, vmax=0.3)
	normalize2 = colors.Normalize(vmin=0, vmax=0.1)
	normalize3 = colors.Normalize(vmin=0, vmax=0.05)
#	plt.scatter(x,y,z,s=5, cmap=colormap, norm=normalize, marker='*')

	fig = plt.figure(figsize=(15,6))
	ax1 = fig.add_subplot(231)
	ax2 = fig.add_subplot(232)
	ax3 = fig.add_subplot(233)
	
	ax4 = fig.add_subplot(234)
	ax5 = fig.add_subplot(235)
	ax6 = fig.add_subplot(236)
	
	def reduce_std(x):
		return np.median(np.abs(x-np.median(x)))

	# Loop through the CBV areas:
	# - or run them in parallel - whatever you like!
	for ii, cbv_area in enumerate(cbv_areas):	
		
#		if ii>5:
#			continue
		
		
		results = np.load('mat-sector%02d-%d_free_weights.npz' % (sector, cbv_area))['res']
		pos_mag = np.zeros([results.shape[0], 3])
		
		for jj, star in enumerate(results[:,0]):
			cursor.execute("""SELECT todolist.starid,todolist.priority,todolist.tmag,eclon,eclat,mean_flux,variance FROM todolist LEFT JOIN diagnostics ON todolist.priority=diagnostics.priority WHERE
				datasource='ffi'
				AND status=1
				AND cbv_area=?
				AND todolist.starid=?
			ORDER BY variability ASC LIMIT 30000;""", (cbv_area, int(star)))		
			star_single = cursor.fetchall()
			
			if cbv_area==122:
				if star_single[0]['eclat']>15:
					pos_mag[jj, 0] =  50
					pos_mag[jj, 1] = 10
					pos_mag[jj, 2] = 10
					continue
			pos_mag[jj, 0] = star_single[0]['eclon']
			pos_mag[jj, 1] = star_single[0]['eclat']
			pos_mag[jj, 2] = star_single[0]['tmag']
#			print(star_single[0]['starid'])
#			print(star_single[0]['eclon'])
#			print(star_single[0]['eclat'])
#			print(star_single[0]['tmag'])
#		
#			if jj==10:
#				sys.exit()

		results=np.column_stack((results, pos_mag))
		print(results.shape, np.min(results[:,1]), np.max(results[:,1]), np.median(np.abs(results[:,1])), np.median(np.abs(results[:,2])), np.median(np.abs(results[:,3])))
		

		ax1.hexbin(results[:,-3], results[:,-1], C=np.abs(results[:,1]), gridsize=10, reduce_C_function=np.median, cmap=colormap, norm=normalize1, marginals=False)
		ax2.hexbin(results[:,-3], results[:,-2], C=np.abs(results[:,2]), gridsize=10, reduce_C_function=np.median, cmap=colormap, norm=normalize2, marginals=False)
		ax3.hexbin(results[:,-3], results[:,-2], C=np.abs(results[:,3]), gridsize=10, reduce_C_function=np.median, cmap=colormap, norm=normalize3, marginals=False)
		
		ax4.hexbin(results[:,-3], results[:,-2], C=np.abs(results[:,1]), gridsize=10, reduce_C_function=reduce_std, cmap=colormap, norm=normalize1, marginals=False)
		ax5.hexbin(results[:,-3], results[:,-2], C=np.abs(results[:,2]), gridsize=10, reduce_C_function=reduce_std, cmap=colormap, norm=normalize2, marginals=False)
		ax6.hexbin(results[:,-3], results[:,-2], C=np.abs(results[:,3]), gridsize=10, reduce_C_function=reduce_std, cmap=colormap, norm=normalize3, marginals=False)

	plt.show()
















