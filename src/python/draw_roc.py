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
from collections import OrderedDict

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
from PIL import Image
from io import BytesIO
import matplotlib.font_manager as font_manager
from matplotlib.pyplot import gca

merged_stats = list()
def MergeStats(work_dirs, grps, grp_prefixes, vmin, vmax, germ, som, vartypes=None):
	global merged_stats
	if som:
		if re.search(',',work_dirs) and re.search('/',grp_prefixes):
			for grp,work_dir,prefixes in zip(grps.split(','),work_dirs.split(','),grp_prefixes.split('/')):
				merged_stats.append([grp,StatsByGrpSomatic(work_dir,prefixes, vmin, vmax)])
		elif not re.search(',',work_dirs) and not re.search('/',grp_prefixes):
			merged_stats.append([grps,StatsByGrpSomatic(work_dirs,grp_prefixes, vmin, vmax)])
		else:
			print "[ERROR] --prefix: number of prefixes shoule match the number of plots!"
			
	elif germ:
		if re.search(',',work_dirs) and re.search('/',grp_prefixes):
			for grp,work_dir,prefixes,vartype in zip(grps.split(','),work_dirs.split(','),grp_prefixes.split('/'),vartypes.split(',')):
				merged_stats.append([grp,StatsByGrpGermline(work_dir,prefixes, vmin, vmax,vartype)])
		elif not re.search(',',work_dirs) and not re.search('/',grp_prefixes):
			merged_stats.append([grps,StatsByGrpGermline(work_dirs,grp_prefixes, vmin, vmax,vartypes)])
		else:
			print "[ERROR] --prefix: number of prefixes shoule match the number of plots!"
			
	else:
		print("Something went wrong with --somatic/--germline]")
		exit(1)
	

def StatsByGrpSomatic(work_dir,prefixes, vmin, vmax):
	global merged_stats
	dic_df_stats = OrderedDict()
	if re.search(',',prefixes):
		for pre_prefix in prefixes.split(','):
			if re.search(':',pre_prefix):
				prefix = pre_prefix.split(':')
			else:
				prefix = [pre_prefix,pre_prefix]
			l_stats = list()
			for filename in glob.glob(work_dir+'/'+prefix[0]+"*.stats.csv"):
				score = float(re.search('[0-9.]+',re.sub("\.stats\.csv","",filename.split('/')[-1])).group().strip('.'))
				if vmin<=score<=vmax:
					df_csv = pd.read_csv(filename)
					l_stats.append([score]+list(df_csv.loc[0,:][['precision','recall']]))
					
			dic_df_stats[prefix[1]] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])
			
	else:
		if re.search(':',prefixes):
			prefix = prefixes.split(':')
		else:
			prefix = [prefixes,prefixes]
		l_stats = list()
		for filename in glob.glob(work_dir+'/'+prefix[0]+"*.stats.csv"):
			score = float(re.search('[0-9.]+',re.sub("\.stats\.csv","",filename.split('/')[-1])).group().strip('.'))
			if vmin<=score<=vmax:
				df_csv = pd.read_csv(filename)
				l_stats.append([score]+list(df_csv.loc[0,:][['precision','recall']]))
				
		dic_df_stats[prefix[1]] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])
		
	return dic_df_stats

def StatsByGrpGermline(work_dir,prefixes, vmin, vmax,vartype):
	global merged_stats
	dic_df_stats = OrderedDict()
	if re.search(',',prefixes):
		for pre_prefix in prefixes.split(','):
			if re.search(':',pre_prefix):
				prefix = pre_prefix.split(':')
			else:
				prefix = [pre_prefix,pre_prefix]
			
			l_stats = list()
			for filename in glob.glob(work_dir+'/'+prefix[0]+"*.summary.csv"):
				score = float(re.search('[0-9.]+',re.sub("\.summary\.csv","",filename.split('/')[-1])).group().strip('.'))
				if vmin<=score<=vmax:
					df_csv = pd.read_csv(filename)
					if vartype in ["Indel","indel","INDEL"]:
						l_stats.append([score]+list(df_csv.loc[1,:][['METRIC.Precision','METRIC.Recall']]))
					elif vartype in ["snp","snv","Snv","Snp","SNV","SNP"]:
						l_stats.append([score]+list(df_csv.loc[3,:][['METRIC.Precision','METRIC.Recall']]))
					
			dic_df_stats[prefix[1]] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])
			
	else:
		if re.search(':',prefixes):
			prefix = prefixes.split(':')
		else:
			prefix = [prefixes,prefixes]
		l_stats = list()
		for filename in glob.glob(work_dir+'/'+prefix[0]+"*.summary.csv"):
			score = float(re.search('[0-9.]+',re.sub("\.summary\.csv","",filename.split('/')[-1])).group().strip('.'))
			if vmin<=score<=vmax:
				df_csv = pd.read_csv(filename)
				if vartype in ["Indel","indel","INDEL"]:
					l_stats.append([score]+list(df_csv.loc[1,:][['METRIC.Precision','METRIC.Recall']]))
				elif vartype in ["snp","snv","Snv","Snp","SNV","SNP"]:
					l_stats.append([score]+list(df_csv.loc[3,:][['METRIC.Precision','METRIC.Recall']]))
				
		dic_df_stats[prefix[1]] = pd.DataFrame(l_stats,columns=['score','precision','recall']).sort_values(by=['score']).set_index(['score'])

	return dic_df_stats
	
	
def DrawROC(output,cols,xmin,xmax,ymin,ymax,num_prefixes,isGray,lcolr,ltype,axisvs,legendloc,font,linewidth,markersize,
legendfontsize,titlefontsize,tickfontsize,labelfontsize,
pre_newxlabels,pre_newylabels):
	global merged_stats
	if isGray:
		plt.style.use('grayscale')

	if re.search(',',legendloc):
		lg_locs = [int(val) for val in legendloc.split(',')]
	elif legendloc == "":
		lg_locs = [0 for i in range(len(merged_stats))]
	else:
		lg_locs = [int(legendloc)]

	if tickfontsize==None:
		tickfont = {'fontname':'serif'}
	else:
		tickfont = {'fontname':'serif','fontsize':tickfontsize}
	if font != "" and legendfontsize!="":
		hfont = {'fontname':font}
		afont = {'fontname':font}
		lfont = font_manager.FontProperties(family=font, weight='normal',style='normal',size=int(legendfontsize))
	elif font != "":
		hfont = {'fontname':font}
		afont = {'fontname':font}
		lfont = font_manager.FontProperties(family=font, weight='normal',style='normal')
	elif legendfontsize!="":
		hfont = {'fontname':'serif'}
		afont = {'fontname':'serif'}
		lfont = font_manager.FontProperties(family='serif',weight='normal',style='normal',size=int(legendfontsize))
	else:
		hfont = {'fontname':'serif'}
		afont = {'fontname':'serif'}
		lfont = font_manager.FontProperties(family='serif',weight='normal',style='normal')

	# border axis visibility
	axisvss = list()
	if re.search(',',axisvs):
		for val in axisvs.split(','):
			if int(val) in [0,1]:
				axisvss.append(int(val))
			else:
				print "input proper value in axis! 0(false) or 1(true) (left,bottom,right,top)"
				exit(1)
		if len(axisvss)!=4:
			print "input proper value in axis! (left,bottom,right,top)"
			exit(1)
	elif axisvs=="":
		axisvss = [1,1,1,1]
	else:
		print "input proper value in axis! (left,bottom,right,top)"
		exit(1)

	# line color
	if re.search(',',lcolr):
		lcolrs = [val for val in lcolr.split(',')]
	elif lcolr == "":
		lcolrs = [None for i in range(num_prefixes)]
	else:
		lcolrs = [lcolr]

	# line stlye
	if re.search(',',ltype):
		ltypes = [val for val in ltype.split(',')]
	elif lcolr == "":
		ltypes = [None for i in range(num_prefixes)]
	else:
		ltypes = [ltype]

	# xlim ylim
	if re.search(',',xmin):
		xmins = [float(num) for num in xmin.split(',')]
	elif xmin == "":
		xmins = [None for i in range(len(merged_stats))]
	else:
		xmins = [float(xmin)]
	if re.search(',',xmax):
		xmaxs = [float(num) for num in xmax.split(',')]
	elif xmax == "":
		xmaxs = [None for i in range(len(merged_stats))]
	else:
		xmaxs = [float(xmax)]
	if re.search(',',ymin):
		ymins = [float(num) for num in ymin.split(',')]
	elif ymin == "":
		ymins = [None for i in range(len(merged_stats))]
	else:
		ymins = [float(ymin)]
	if re.search(',',ymax):
		ymaxs = [float(num) for num in ymax.split(',')]
	elif ymax == "":
		ymaxs = [None for i in range(len(merged_stats))]
	else:
		ymaxs = [float(ymax)]

	if pre_newxlabels in [None,""]:
		newxlabels = ["" for i in range(len(merged_stats))]
	elif re.search(',',pre_newxlabels):
		newxlabels = pre_newxlabels.split(',')
	else:
		print "label should be splited by ','!"
	if pre_newylabels in [None,""]:
		newylabels = ["" for i in range(len(merged_stats))]
	elif re.search(',',pre_newylabels):
		newylabels = pre_newylabels.split(',')
	else:
		print "label should be splited by ','!"

	total = 0
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
		for [grp,stats],xmin,xmax,ymin,ymax,lg_loc,newxlabel,newylabel in zip(merged_stats,xmins,xmaxs,ymins,ymaxs,lg_locs,newxlabels,newylabels):
			l_lgn = list()
			for [key,stat],lcolr,ltype in zip(stats.items(),lcolrs,ltypes):
				axes.plot(stat['recall'],stat['precision'],c=lcolr,ls=ltype,lw=float(linewidth))
				axes.scatter(stat['recall'].iloc[0],stat['precision'].iloc[0],c=lcolr,s=markersize)
				l_lgn.append(key)
			axes.legend(l_lgn,prop=lfont,loc=lg_loc)
			if titlefontsize!=None:
				axes.set_title(grp,fontsize=titlefontsize,**hfont)
			else:
				axes.set_title(grp,**hfont)

			if xmin==None:
				xmin=stat['recall'].iloc[0]-0.3
			if xmax==None:
				xmax=stat['recall'].iloc[0]+0.1
			if ymin==None:
				ymin=stat['precision'].iloc[0]-0.2
			if ymax==None:
				ymax=stat['precision'].iloc[0]+0.2

			if labelfontsize!=None:
				if newxlabel == "":
					axes.xaxis.label.set_visible(False)
				else:
					axes.set_xlabel(newxlabel,fontsize=labelfontsize,**afont)
				if newylabel == "":
					axes.yaxis.label.set_visible(False)
				else:
					axes.set_ylabel(newylabel,fontsize=labelfontsize,**afont)
			else:
				if newxlabel == "":
					axes.xaxis.label.set_visible(False)
				else:
					axes.set_xlabel(newxlabel,**afont)
				if newylabel == "":
					axes.yaxis.label.set_visible(False)
				else:
					axes.set_ylabel(newylabel,**afont)

			# axes.set_xticks(axes.get_xticks())
			# axes.set_yticks(axes.get_yticks())
			# if tickfontsize!=None:
			# 	axes.set_xticklabels(axes.get_xticks(),tickfont)
			# 	axes.set_yticklabels(axes.get_yticks(),tickfont)
			# else:
			# 	axes.set_xticklabels(axes.get_xticks(), afont)
			# 	axes.set_yticklabels(axes.get_yticks(), afont)
			axes.set_xticks([xmin,xmin+(xmax-xmin)/2,xmax])
			axes.set_yticks([ymin,ymin+(ymax-ymin)/2,ymax])
			if tickfontsize!=None:
				axes.set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax],tickfont)
				axes.set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax],tickfont)
			else:
				axes.set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax], afont)
				axes.set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax], afont)

			axes.set_xlim([xmin,xmax])
			axes.set_ylim([ymin,ymax])

			axes.spines['left'].set_visible(axisvss[0])
			axes.spines['bottom'].set_visible(axisvss[1])
			axes.spines['right'].set_visible(axisvss[2])
			axes.spines['top'].set_visible(axisvss[3])

	elif rows == 1:
		plt.rcParams["font.family"] = "serif"
		for [grp,stats],xmin,xmax,ymin,ymax,lg_loc,newxlabel,newylabel in zip(merged_stats,xmins,xmaxs,ymins,ymaxs,lg_locs,newxlabels,newylabels):
			l_lgn = list()
			for [key,stat],lcolr,ltype in zip(stats.items(),lcolrs,ltypes):
				axes[col].plot(stat['recall'],stat['precision'],c=lcolr,ls=ltype,lw=float(linewidth))
				axes[col].scatter(stat['recall'].iloc[0],stat['precision'].iloc[0],c=lcolr,s=markersize)
				l_lgn.append(key)
			axes[col].legend(l_lgn,prop=lfont,loc=lg_loc)
			if titlefontsize!=None:
				axes[col].set_title(grp,fontsize=titlefontsize,**hfont)
			else:
				axes[col].set_title(grp,**hfont)

			if xmin==None:
				xmin=stat['recall'].iloc[0]-0.3
			if xmax==None:
				xmax=stat['recall'].iloc[0]+0.1
			if ymin==None:
				ymin=stat['precision'].iloc[0]-0.2
			if ymax==None:
				ymax=stat['precision'].iloc[0]+0.2

			
			if labelfontsize!=None:
				if newxlabel == "":
					axes[col].xaxis.label.set_visible(False)
				else:
					axes[col].set_xlabel(newxlabel,fontsize=labelfontsize,**afont)
				if newylabel == "":
					axes[col].yaxis.label.set_visible(False)
				else:
					axes[col].set_ylabel(newylabel,fontsize=labelfontsize,**afont)
			else:
				if newxlabel == "":
					axes[col].xaxis.label.set_visible(False)
				else:
					axes[col].set_xlabel(newxlabel,**afont)
				if newylabel == "":
					axes[col].yaxis.label.set_visible(False)
				else:
					axes[col].set_ylabel(newylabel,**afont)				

			
			# axes[col].set_xticks(axes[col].get_xticks())
			# axes[col].set_yticks(axes[col].get_yticks())
			# if tickfontsize!=None:
			# 	axes[col].set_xticklabels(axes[col].get_xticks(),tickfont)
			# 	axes[col].set_yticklabels(axes[col].get_yticks(),tickfont)
			# else:
			# 	axes[col].set_xticklabels(axes[col].get_xticks(), afont)
			# 	axes[col].set_yticklabels(axes[col].get_yticks(), afont)

			axes[col].set_xticks([xmin,xmin+(xmax-xmin)/2,xmax])
			axes[col].set_yticks([ymin,ymin+(ymax-ymin)/2,ymax])
			if tickfontsize!=None:
				axes[col].set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax],tickfont)
				axes[col].set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax],tickfont)
			else:
				axes[col].set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax], afont)
				axes[col].set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax], afont)

			axes[col].set_xlim([xmin,xmax])
			axes[col].set_ylim([ymin,ymax])

			axes[col].spines['left'].set_visible(axisvss[0])
			axes[col].spines['bottom'].set_visible(axisvss[1])
			axes[col].spines['right'].set_visible(axisvss[2])
			axes[col].spines['top'].set_visible(axisvss[3])

			col += 1
			
	else:
		for [grp,stats],xmin,xmax,ymin,ymax,lg_loc,newxlabel,newylabel in zip(merged_stats,xmins,xmaxs,ymins,ymaxs,lg_locs,newxlabels,newylabels):
			l_lgn = list()
			for [key,stat],lcolr,ltype in zip(stats.items(),lcolrs,ltypes):
				axes[row,col].plot(stat['recall'],stat['precision'],c=lcolr,ls=ltype,lw=float(linewidth))
				axes[row,col].scatter(stat['recall'].iloc[0],stat['precision'].iloc[0],c=lcolr,s=markersize)
				l_lgn.append(key)
			axes[row,col].legend(l_lgn,prop=lfont,loc=lg_loc)
			if titlefontsize!=None:
				axes[row,col].set_title(grp,fontsize=titlefontsize,**hfont)
			else:
				axes[row,col].set_title(grp,**hfont)

			if xmin==None:
				xmin=stat['recall'].iloc[0]-0.3
			if xmax==None:
				xmax=stat['recall'].iloc[0]+0.1
			if ymin==None:
				ymin=stat['precision'].iloc[0]-0.2
			if ymax==None:
				ymax=stat['precision'].iloc[0]+0.2

			if labelfontsize!=None:
				if newxlabel == "":
					axes[row,col].xaxis.label.set_visible(False)
				else:
					axes[row,col].set_xlabel(newxlabel,fontsize=labelfontsize,**afont)
				if newylabel == "":
					axes[row,col].yaxis.label.set_visible(False)
				else:
					axes[row,col].set_ylabel(newylabel,fontsize=labelfontsize,**afont)
			else:
				if newxlabel == "":
					axes[row,col].xaxis.label.set_visible(False)
				else:
					axes[row,col].set_xlabel(newxlabel,**afont)
				if newylabel == "":
					axes[row,col].yaxis.label.set_visible(False)
				else:
					axes[row,col].set_ylabel(newylabel,**afont)
				
			# axes[row,col].set_xticks(axes[row,col].get_xticks())
			# axes[row,col].set_yticks(axes[row,col].get_yticks())
			# if tickfontsize!=None:
			# 	axes[row,col].set_xticklabels(axes[row,col].get_xticks(),tickfont)
			# 	axes[row,col].set_yticklabels(axes[row,col].get_yticks(),tickfont)
			# else:
			# 	axes[row,col].set_xticklabels(axes[row,col].get_xticks(), afont)
			# 	axes[row,col].set_yticklabels(axes[row,col].get_yticks(), afont)
			axes[row,col].set_xticks([xmin,xmin+(xmax-xmin)/2,xmax])
			axes[row,col].set_yticks([ymin,ymin+(ymax-ymin)/2,ymax])
			if tickfontsize!=None:
				axes[row,col].set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax],tickfont)
				axes[row,col].set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax],tickfont)
			else:
				axes[row,col].set_xticklabels([xmin,xmin+(xmax-xmin)/2,xmax], afont)
				axes[row,col].set_yticklabels([ymin,ymin+(ymax-ymin)/2,ymax], afont)

			axes[row,col].set_xlim([xmin,xmax])
			axes[row,col].set_ylim([ymin,ymax])

			axes[row,col].spines['left'].set_visible(axisvss[0])
			axes[row,col].spines['bottom'].set_visible(axisvss[1])
			axes[row,col].spines['right'].set_visible(axisvss[2])
			axes[row,col].spines['top'].set_visible(axisvss[3])

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
	parser.add_argument('--type',dest="type",type=str,help='indel/snv',default=None,required=False)
	parser.add_argument('--out',type=str,help=' ',required=True)
	parser.add_argument('--col',dest="cols",type=int,help=' ',required=True)
	parser.add_argument('--min',dest="vmin",type=float,help=' ',default=-999999,required=False)
	parser.add_argument('--max',dest="vmax",type=float,help=' ',default=999999,required=False)
	parser.add_argument('--germline',dest="germ",action='store_true')
	parser.add_argument('--somatic',dest="som",action='store_true')
	parser.add_argument('--gray',dest="gray",action='store_true')
	parser.add_argument('--xmin',dest="xmin",type=str,help=' ',default="",required=False)
	parser.add_argument('--xmax',dest="xmax",type=str,help=' ',default="",required=False)
	parser.add_argument('--ymin',dest="ymin",type=str,help=' ',default="",required=False)
	parser.add_argument('--ymax',dest="ymax",type=str,help=' ',default="",required=False)

	parser.add_argument('--lcolr',dest="lcolr",type=str,help=' ',default="",required=False)
	parser.add_argument('--ltype',dest="ltype",type=str,help=' ',default="",required=False)
	parser.add_argument('--axis',dest="axisvs",type=str,help=' ',default="",required=False)
	parser.add_argument('--legendloc',dest="legendloc",type=str,help=' ',default="",required=False)
	parser.add_argument('--font',dest="font",type=str,help=' ',default="serif",required=False)
	parser.add_argument('--linewidth',dest="linewidth",type=float,help=' ',default=1,required=False)
	parser.add_argument('--markersize',dest="markersize",type=float,help=' ',default=5,required=False)

	parser.add_argument('--legendfontsize',dest="legendfontsize",type=float,help=' ',default=None,required=False)
	parser.add_argument('--titlefontsize',dest="titlefontsize",type=float,help=' ',default=None,required=False)
	parser.add_argument('--tickfontsize',dest="tickfontsize",type=float,help=' ',default=None,required=False)
	parser.add_argument('--labelfontsize',dest="labelfontsize",type=float,help=' ',default=None,required=False)

	parser.add_argument('--xlabel',dest="newxlabel",type=str,help=' ',default=None,required=False)
	parser.add_argument('--ylabel',dest="newylabel",type=str,help=' ',default=None,required=False)
	
	args = parser.parse_args()
	try:
		if not args.germ and not args.som:
			print("Either of --somatic or --germline should be used!")
			exit(1)
		if args.germ and args.som:
			print("Only one can be chosen between --somatic or --germline!")
			exit(1)
		if args.germ and args.type==None:
			print("--type is needed for plotting the germline!")
			exit(1)
		MergeStats(args.dir,args.grp,args.prefix, args.vmin, args.vmax, args.germ, args.som,args.type)

		num_prefixes = len(args.prefix.split(','))
		DrawROC(args.out,args.cols,args.xmin,args.xmax,args.ymin,args.ymax,num_prefixes,
		args.gray,args.lcolr,args.ltype,args.axisvs,args.legendloc,args.font,args.linewidth,args.markersize,
		args.legendfontsize,args.titlefontsize,args.tickfontsize,args.labelfontsize,
		args.newxlabel,args.newylabel)
			
	except Exception as e:
		logger.error(traceback.format_exc())
		raise e