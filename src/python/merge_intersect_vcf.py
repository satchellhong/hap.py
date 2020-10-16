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

class MergeVcf:
	def __init__(self):
		self.dict_merged = dict()
		self.keysets = list()
		self.intersect = list()

	def MergeVcfs(self,vcfs, out,isSomatic=False):
		isHg38 = False
		meta_data = list()

		if re.search(',',vcfs):
			isFirst = True
			for vcf in vcfs.split(','):
				if isFirst:
					with open(vcf,'r') as vcf_in:
						for line in vcf_in:
							if line.startswith('#'):
								meta_data.append(line.strip())
							else:
								if line.startswith('chr'):
									isHg38 = True
								break
				self.ReadVcf(vcf,isFirst)
				isFirst = False
		else:
			print "Only one vcf file cannot be merged with itself.\nAt least two files are needed."

		self.intersect = set.intersection(*self.keysets)

		print "Organizing merged VCF data..."

		meta_data[-1] = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT"

		df_merged = DataFrame(self.dict_merged.values()).astype({0:"category",1:int})
		df_merged.set_index(9,inplace=True)
		df_merged = df_merged.loc[self.intersect,:].sort_values(by=[0,1])
		df_merged = df_merged.astype({0:str,1:str})

		if isHg38:
			df_merged[0] = "chr"+df_merged[0]

		df_values = df_merged.iloc[:,:9].values
		
		print "Writing to "+out
		with open(out,'w') as vcf_out:
			for line in meta_data:
				vcf_out.write(line+'\n')
			for line in df_values:
				vcf_out.write('\t'.join(line)+'\n')

	def ReadVcf(self,vcf,isFirst):
		print "[Merging VCF]: " + vcf

		with open(vcf,'r') as vcf_in:
			temp = list()
			for line in vcf_in:
				if not line.startswith('#'):
					elements = line.strip().split()
					sKey = elements[0]+'_'+elements[1]
					temp.append(sKey)

					if self.dict_merged.get(sKey) == None:
						self.dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:9]+[sKey]
					else:
						sRef1 = self.dict_merged[sKey][3]
						sRef2 = elements[3]
						sAlt1 = self.dict_merged[sKey][4]
						sAlt2 = elements[4]

						if sRef1 == sRef2 and sAlt1 == sAlt2:
							continue

						if sRef1 != sRef2:
							if re.search(',',sRef1):
								lRef1 = set(sRef1.split(','))
							else:
								lRef1 = set(sRef1)
							if re.search(',',sRef2):
								lRef2 = set(sRef2.split(','))
							else:
								lRef2 = set(sRef2)
							lRef = lRef1.intersection(lRef2)
							if len(lRef)>0:
								elements[3] = ','.join(lRef)
							else:
								continue

						if sAlt1 != sAlt2:
							if re.search(',',sAlt1):
								lAlt1 = set(sAlt1.split(','))
							else:
								lAlt1 = set(sAlt1)
							if re.search(',',sAlt2):
								lAlt2 = set(sAlt2.split(','))
							else:
								lAlt2 = set(sAlt2)
							lAlt = lAlt1.intersection(lAlt2)
							if len(lAlt)>0:
								elements[4] = ','.join(lAlt)
							else:
								continue

						self.dict_merged[sKey] = [re.sub("chr","",elements[0])]+elements[1:9]+[sKey]

			self.keysets.append(set(temp))


if __name__ == "__main__":  
	FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
	logging.basicConfig(level=logging.INFO, format=FORMAT)
	logger = logging.getLogger(__name__)
	
	parser = argparse.ArgumentParser(description='Merge vcf files')
	
	parser.add_argument('--vcfs',type=str,help=' ',required=True)
	parser.add_argument('--out',dest="out",type=str,help='output filename',default=None,required=True)
	
	args = parser.parse_args()
	
	try:
		MergeVcf = MergeVcf()
		MergeVcf.MergeVcfs(args.vcfs,args.out)
			
	except Exception as e:
		logger.error(traceback.format_exc())
		raise e