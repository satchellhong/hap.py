#!/usr/bin/env python
# coding=utf-8

import os
import time
import argparse
import logging
import traceback
import pandas as pd
import numpy as np
from pandas import DataFrame
from numpy.random import *
import glob
import re

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from PIL import Image
from io import BytesIO

merged_stats = list()
def MergeStats(work_dirs, grps, prefixes, vmin, vmax):
	if re.search(',',work_dirs):
		for grp,work_dir in zip(grps.split(','),work_dirs.split(',')):
			merged_stats.append([grp,StatsByGrp(work_dir,prefixes, vmin, vmax)])
	else:
		merged_stats.append([grps,StatsByGrp(work_dirs,prefixes, vmin, vmax)])

def StatsByGrp(work_dir,prefixes, vmin, vmax):
	dic_df_stats = dict()
	if re.search(',',prefixes):
		for prefix in prefixes.split(','):
			l_stats = list()
			for filename in glob.glob(work_dir+'/'+prefix+"*.stats.csv"):
				score = float(re.search('[0-9.]+',re.sub("\.stats\.csv","",filename.split('/')[-1])).group().strip('.'))
				if vmin<=score<=vmax:
					df_csv = pd.read_csv(filename)
					l_stats.append([score]+list(df_csv.loc[0,:][['precision','recall']]))
					
			dic_df_stats[prefix] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])
			
	else:
		l_stats = list()
		for filename in glob.glob(work_dir+'/'+prefixes+"*.stats.csv"):
			score = float(re.search('[0-9.]+',re.sub("\.stats\.csv","",filename.split('/')[-1])).group().strip('.'))
			if vmin<=score<=vmax:
				df_csv = pd.read_csv(filename)
				l_stats.append([score]+list(df_csv.loc[0,:][['precision','recall']]))
				
		dic_df_stats[prefixes] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])
		
	return dic_df_stats
	
	
def DrawROC(output,cols):
	
	total = 0;
	for grp,stats in merged_stats:
		print grp
		for key,stat in stats.items():
			print key
			print stat
		
	rows = 1
	if cols>=len(merged_stats):
		cols = len(merged_stats)
	else:
		rows = (len(merged_stats)/cols)
		if len(merged_stats)%cols!=0:
			rows+=1
	
	figure, axes = plt.subplots(nrows=rows, ncols=cols, figsize = (cols*3,rows*3))
	
	row = 0
	col = 0
	if rows == 1 and col == 1:
		for grp,stats in merged_stats:
			l_lgn = list()
			for key,stat in stats.items():
				axes.plot(stat['recall'],stat['precision'])
				axes.scatter(stat.loc[0,'recall'],stat.loc[0,'precision'])
				l_lgn.append(key)
			axes.legend(l_lgn)
			axes.set_title(grp)
	elif rows == 1:
		for grp,stats in merged_stats:
			l_lgn = list()
			for key,stat in stats.items():
				axes[col].plot(stat['recall'],stat['precision'])
				axes[col].scatter(stat.loc[0,'recall'],stat.loc[0,'precision'])
				l_lgn.append(key)
			axes[col].legend(l_lgn)
			axes[col].set_title(grp)
			col += 1
	else:
		for grp,stats in merged_stats:
			l_lgn = list()
			for key,stat in stats.items():
				axes[row,col].plot(stat['recall'],stat['precision'])
				axes[row,col].scatter(stat.loc[0,'recall'],stat.loc[0,'precision'])
				l_lgn.append(key)
			axes[row,col].legend(l_lgn)
			axes[row,col].set_title(grp)
			if col<cols:
				col += 1
			else:
				col = 0
				row += 1
		
	gs = gridspec.GridSpec(1, 1)
	figure.tight_layout(w_pad=0.5, h_pad=1.0, rect=[0,0,1,0.985])
	
	png1 = BytesIO()
	plt.savefig(output+'.png',format='png')
	
	png2 = Image.open(output+'.png')
	png2.save(output+'.tiff')
	png1.close()
	
if __name__ == "__main__":  
	FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
	logging.basicConfig(level=logging.INFO, format=FORMAT)
	logger = logging.getLogger(__name__)
	
	parser = argparse.ArgumentParser(description='Draw PR curve')
	
	parser.add_argument('--dir',type=str,help=' ',required=True)
	parser.add_argument('--grp',type=str,help=' ',required=True)
	parser.add_argument('--prefix',type=str,help=' ',required=True)
	parser.add_argument('--out',type=str,help=' ',required=True)
	parser.add_argument('--col',dest="cols",type=int,help=' ',required=True)
	parser.add_argument('--min',dest="vmin",type=float,help=' ',default=-999999,required=False)
	parser.add_argument('--max',dest="vmax",type=float,help=' ',default=999999,required=False)
	
	
	args = parser.parse_args()
	
	try:
		MergeStats(args.dir,args.grp,args.prefix, args.vmin, args.vmax)
		DrawROC(args.out,args.cols)
			
	except Exception as e:
		logger.error(traceback.format_exc())
		raise e