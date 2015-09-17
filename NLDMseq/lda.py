#!/usr/bin/env python
#coding=utf-8
import os
import os.path
import shutil
import glob
import ctypes
import multiprocessing as mp
########################################
#function:copy file from one dir to another
#input: res file path, target file path
#output:none
########################################
def copy_file(res_file,tar_file):
	shutil.copy(res_file,tar_file)
########################################
#function:remove all file in the dir and the floder itself
#input: dir to remove
#output:none
########################################
def removeFileInFirstDir(targetDir): 
	for file_path in os.listdir(targetDir): 
		targetFile = os.path.join(targetDir,file_path) 
		if os.path.isfile(targetFile): 
			os.remove(targetFile)
        		
########################################
#function: split the gene and isoform abundance according the annotation
#input: refFlat_Check_GeneUnionExonLen ,expression_log
#output: gene_expre isoform_expre
########################################
def gene_expre_log(data_path,output_path):
	#check the output file exists
	gene_output_file_path =os.path.join(output_path,'expression_gene.txt')
	iso_output_file_path  =os.path.join(output_path,'expression_isoform.txt')

	if(os.path.exists(gene_output_file_path)):
		os.remove(gene_output_file_path)
	if(os.path.exists(iso_output_file_path)):
		os.remove(iso_output_file_path)

	annotation=os.path.join(data_path,'NLDMseq.Gene.Info')
	log_path  =os.path.join(data_path,'expression_log')
	ann={}
	#read annotation txet and save gene:isoform to a dictionary
	fp=open(annotation,'r')
	for line in fp:
		line=line.rstrip()
		line=line.split('\t')
		ann[line[0]]=line[3:len(line)]
	fp.close()
        #split the expression_log to gene_expre and isoform_expre
	
	gene_log=open(gene_output_file_path,'w')
	isoform_log=open(iso_output_file_path,'w')
	gene_file_list=os.listdir(log_path);
	for gene_file in gene_file_list:
		fp=open(os.path.join(log_path,gene_file),'r')
		for line in fp:
			line=line.rstrip()
			line=line.split('\t')
			gene_log.write(line[0]+'\t')
			gene_log.write(line[len(line)-2]+'\t')
			gene_log.write(line[len(line)-3]+'\n')
			t_cnt=0
			for item in ann[line[0]]:
				isoform_log.write(item+'\t'+line[0]+'\t'+line[t_cnt+1]+'\t'+line[t_cnt+2]+'\n')
				t_cnt=t_cnt+1
	gene_log.close()
	isoform_log.close()
########################################
#function: calcultate the abundances of gene and isoform 
#input:EM parement
#output:count of gene
########################################
def pp_lda(gene_list,isoLen_path,Map_path,NormData_path,expression_log_path,data_path,SEQ_DEPTH):
	for item in gene_list:
		doc_str=item
		isoLen_exists   =os.path.exists(os.path.join(isoLen_path,doc_str))
		map_exists      =os.path.exists(os.path.join(Map_path,doc_str))
		normdata_exists =os.path.exists(os.path.join(NormData_path,doc_str))
		if(isoLen_exists and map_exists and normdata_exists):
			expression_log_out=os.path.join(expression_log_path,item)
			so.gene_expre(data_path,doc_str,SEQ_DEPTH,expression_log_out)

def gene_expre(data_path,SEQ_DEPTH):
	isoLen_path  =  os.path.join(data_path,'isoLen')
	Data_path    =  os.path.join(data_path,'ModelMultiGene_Data')
	Map_path     =  os.path.join(data_path,'ModelMultiGene_Map')
	NormData_path=  os.path.join(data_path,'ModelMultiGene_NormData')
	expression_log_path=os.path.join(data_path,'expression_log')
	if(os.path.exists(expression_log_path)):
		shutil.rmtree(expression_log_path)
	os.mkdir(expression_log_path)
	
	res_list=os.listdir(Data_path)

	res_list=os.listdir(Data_path)
	block_size=1000
	start=0
	end=len(res_list)
	endi=0;
	starti=0;
	pool=mp.Pool(processes=mp.cpu_count()*2)
	while(endi<end):
		endt=starti+block_size;
		endi=min(endt,end);
		pool.apply_async(pp_lda,(res_list[starti:endi],isoLen_path,Map_path,NormData_path,expression_log_path,data_path,SEQ_DEPTH, ))
		starti=starti+block_size;
	pool.close()
	pool.join()
########################################
#function:lda main 
#input: data input ,Sequence depth ,expression output path
#output:none
########################################
so=ctypes.CDLL('./libgene_expre.so')
def lda_expre(data_path,SEQ_DEPTH,output_file):
	gene_expre(data_path,SEQ_DEPTH)
	gene_expre_log(data_path,output_file)

        

