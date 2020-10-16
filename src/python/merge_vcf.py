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
from collections import Counter

dict_merged = dict()

def MergeVcfs(vcfs, out,isSomatic=False):
	global dict_merged

	isHg38 = False
	meta_data = list()

	if re.search(',',vcfs):
		isFirst = True
		for vcf in vcfs.split(','):
			with open(vcf,'r') as vcf_in:
				for line in vcf_in:
					if isFirst:
						if line.startswith('#'):
							meta_data.append(line.strip())
						else:
							if line.startswith('chr'):
								isHg38 = True
							break

			ReadVcf(vcf,isSomatic)
			isFirst = False
	else:
		print "Only one vcf file cannot be merged with itself.\nAt least two files are needed."

	if isSomatic:
		meta_data[-1] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"
	else:
		meta_data[-1] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"

	print "Sorting merged VCF data..."

	df_merged = DataFrame(dict_merged.values()).astype({0:"category",1:int}).sort_values(by=[0,1])
	df_merged = df_merged.astype({0:str,1:str})
	
	print "Organizing merged VCF data..."

	if isHg38:
		df_merged[0] = "chr"+df_merged[0]
	
	df_values = None
	if isSomatic:
		df_values = df_merged.iloc[:,:9].values
	else:
		df_values = df_merged.iloc[:,:10].values
		
		# 각 row 별로 어떤 genotype(0/0, 0/1, 1/1, 1/2,...)이 가장 많이 나왔는지 voting해서 가장 많이 나온 것을 선택한다 
		# find  most frequent element in the list
		for row_idx,line in enumerate(df_values):
			if re.search(':',line[9]):
				genos = re.sub('\|','/',line[9]).split(':')
				# print max(set(genos), key = genos.count) # 가장 frequent한 값을 return해줌
				# print Counter(genos).most_common() #각 element frequency를 내림차순으로 보여줌

				selected = None
				for gen_idx,gen in enumerate(genos):
					if gen == max(set(genos), key = genos.count):
						selected=gen_idx
						line[9] = genos[selected]
						break

				line[3] = line[3].split(':')[selected]
				line[4] = line[4].split(':')[selected]
				# df_values[row_idx] = line # no need because of shallow copy

	print "Writing to "+out
	with open(out,'w') as vcf_out:
		for line in meta_data:
			vcf_out.write(line+'\n')
		for line in df_values:
			vcf_out.write('\t'.join(line)+'\n')


def ReadVcf(vcf,isSomatic):
	global dict_merged
	if isSomatic:
		print "[Merging Somatic VCF]: " + vcf
	else:
		print "[Merging Germline VCF]: " + vcf
		
	with open(vcf,'r') as vcf_in:
		for line in vcf_in:
			if not line.startswith('#'):
				elements = line.strip().split()
				sKey = elements[0]+'_'+elements[1]
				sInfo = elements[7]
				fRdscore = 0.0
				if re.search("RDscore=[0-9\.]+",sInfo):
					fRdscore = float(re.search("RDscore=[0-9\.]+",sInfo).group().split('=')[-1])
				else:
					print "[WARNING]: {0} : {1} {2} has no RDscore value!".format(vcf,elements[0],elements[1])
					continue
				
				if dict_merged.get(sKey) == None:
					if isSomatic:
						dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:9]+[fRdscore]
					else:
						sGt = None
						for i,s in enumerate(elements[8].split(':')):
							if s=="GT":
								sGt = elements[9].split(':')[i]
								break
						if sGt == None:
							print "No GT found in format!"
							print vcf
							print elements
							exit(1)

						dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:8]+["GT", sGt, fRdscore]
				else:
					if isSomatic:
						if float(dict_merged[sKey][-1])<fRdscore:
							dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:9]+[fRdscore]
					else:
						if float(dict_merged[sKey][-1])<fRdscore:
							sGt = None
							for i,s in enumerate(elements[8].split(':')):
								if s=="GT":
									sGt = elements[9].split(':')[i]
									break
							if sGt == None:
								print "No GT found in format!"
								print vcf
								print elements
								exit(1)

							dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:8]+["GT", sGt, fRdscore]
						
						if float(dict_merged[sKey][-1])==fRdscore:
							sGt = None
							for i,s in enumerate(elements[8].split(':')):
								if s=="GT":
									sGt = elements[9].split(':')[i]
									break
							if sGt == None:
								print "No GT found in format!"
								print vcf
								print elements
								exit(1)

							elements[3] = dict_merged[sKey][3]+':'+elements[3]
							elements[4] = dict_merged[sKey][4]+':'+elements[4]
							dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:8]+["GT", dict_merged[sKey][9]+':'+sGt, fRdscore]

if __name__ == "__main__":  
	FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
	logging.basicConfig(level=logging.INFO, format=FORMAT)
	logger = logging.getLogger(__name__)
	
	parser = argparse.ArgumentParser(description='Merge vcf files')
	
	parser.add_argument('--vcfs',type=str,help=' ',required=True)
	parser.add_argument('--out',dest="out",type=str,help='output filename',default=None,required=True)
	parser.add_argument('--somatic',dest="som",action='store_true')
	
	args = parser.parse_args()
	
	try:
		MergeVcfs(args.vcfs,args.out,args.som)
			
	except Exception as e:
		logger.error(traceback.format_exc())
		raise e