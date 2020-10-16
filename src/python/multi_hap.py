#!/usr/bin/env python
# coding=utf-8

from subprocess import call, check_output, Popen, PIPE, STDOUT, CalledProcessError
import os
from glob import glob
import argparse
import traceback
import logging
import shutil,psutil
import re
from tqdm import tqdm
from multiprocessing import Pool,cpu_count,freeze_support,RLock
import signal

def Sompy(sompy_cmds,pid=0):
	pbar = tqdm(total=len(sompy_cmds), position=pid+1)
	for sompy_cmd in sompy_cmds:
		pbar.set_description_str("[{0}]".format(sompy_cmd[1]))
		pbar.update(1)
		popen = Popen(sompy_cmd[0], shell=True, stdout=PIPE, universal_newlines=True,stderr=STDOUT)
		popen.communicate()
		if popen.returncode!=0:
			exit(popen.returncode)
	
def MultSompy(truth, query, FP, output, ref, thread=cpu_count()):
	print "Thread: {0}".format(thread)
	
	if not output.endswith('/'):
		output += '/'
	if not os.path.exists(output):
		os.makedirs(output)
			
	vcf_list = list()
	for vcf in glob(query+"*"):
		vcf_list.append(vcf)
		
	# load balancing
	t_vcf_list=list()
	chunk_size = int(len(vcf_list)/thread)
	if chunk_size<1:
		thread=len(vcf_list)
		for vcf in vcf_list:
			out_prefix = output+'/'+vcf.split('/')[-1]
			if(ref==None):
				cmd = "/opt/hap.py/bin/hap.py "+truth+' '+vcf+' -f '+FP+' -o '+out_prefix
			else:
				cmd = "/opt/hap.py/bin/hap.py "+truth+' '+vcf+' -f '+FP+' -r '+ref+' -o '+out_prefix
			t_vcf_list.append([[cmd,vcf.split('/')[-1]]])
			
	else:
		for i,vcf in enumerate(vcf_list):
			out_prefix = output+'/'+vcf.split('/')[-1]
			if(ref==None):
				cmd = "/opt/hap.py/bin/hap.py "+truth+' '+vcf+' -f '+FP+' -o '+out_prefix
			else:
				cmd = "/opt/hap.py/bin/hap.py "+truth+' '+vcf+' -f '+FP+' -r '+ref+' -o '+out_prefix
				
			if i<thread:
				t_vcf_list.append([[cmd,vcf.split('/')[-1]]])
			else:
				t_vcf_list[i%thread].append([cmd,vcf.split('/')[-1]])
	
	
	pool = Pool(processes=thread, initargs=(RLock(),), initializer=init_worker)
	for pid,cmds in enumerate(t_vcf_list):
		pool.apply_async(Sompy, (cmds,pid,))
	
	pool.close()
	pool.join()
		
	print('\n'*(thread+1))

def init_worker(rlock):
	tqdm.set_lock(rlock)
	signal.signal(signal.SIGINT, signal.SIG_IGN)
	
if __name__ == "__main__":
	FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
	logging.basicConfig(level=logging.INFO, format=FORMAT)
	logger = logging.getLogger(__name__)
	
	parser = argparse.ArgumentParser("Somatic Comparison")
	
	parser.add_argument("truth", help="Truth VCF file")
	parser.add_argument("query", help="Query VCF file")
	parser.add_argument("-f", "--false-positives", dest="FP", required=True,
                        help="False-positive region bed file to distinguish UNK from FP")
	parser.add_argument("-o", "--output", dest="output", required=True,
                        help="Output file prefix for statistics and feature table (when selected)")
	parser.add_argument("-r", "--reference", dest="ref", default=None, required=False,
                        help="Specify a reference file.")
	parser.add_argument("--thread", dest="thread",type=int,default=cpu_count(), required=False)
	parser.add_argument("--force", dest="force",action='store_true', default=False)
	
	args = parser.parse_args()
	
	logger.info(args)

	try:
		mem_usage=3
		num_thread = args.thread
		mem_size = psutil.virtual_memory().total/float(1024*1024*1024)
		available_thread = int(mem_size/mem_usage)
		recommended_thread = cpu_count()
		if num_thread>available_thread:
			print("\n\033[93m[[WARNING]]")
			print("RAM size is too small to use that amount of threads : {0}".format(num_thread))
			print("Available number of threads : {0}".format(available_thread))
			num_thread = available_thread
			
			if num_thread>recommended_thread:
				print("Recommended number of threads : {0}".format(recommended_thread))
				if not args.force:
					print("Use --force if you don't want to follow the recommendation.")
					num_thread = recommended_thread
			
			print("Recalibrating the number of threads to {0}\033[0m".format(num_thread))
		
		MultSompy(args.truth,args.query,args.FP,args.output,args.ref,num_thread)
		
	except Exception:
		logger.error(traceback.format_exc())