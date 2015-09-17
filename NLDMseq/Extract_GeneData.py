#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  

import Extract_BowtieData
import static_exonJunc
import static_readOnExon
import geneExonLen
import getIsoLength

import plusGeneName
import Extract_CalcAbsLoc


import os
import os.path

import pp





def startProcess(readFileList,AnnotationDir,ReadLen,readType):
	seq_depth=0
	#Sample = ['Brain','UHR']
	Sample = ['workfloder']
	rep=readFileList
	lanFile=[]
	for i in range(len(readFileList)):
		lanFile.append('Lane'+str(i+1))
	#Create Run path function

	targetDir = AnnotationDir
	workfloder      =  os.path.join(targetDir,'workfloder')

	####Annocation Files Path 
	GeneExonSplit   =  os.path.join(workfloder,'NLDMseq.Exon.Split')
	GeneFile        =  os.path.join(workfloder,'NLDMseq.Gene.Info')
	#####workfloder path
	geneMapout      =  os.path.join(workfloder,'ModelMultiGene_Map')
	geneExonout     =  os.path.join(workfloder,'GeneExonFile')
	GeneData_Path   =  os.path.join(workfloder,'ExtractMultiGeneData')
	GeneAbsLoc_path =  os.path.join(workfloder,'ModelMultiGene_AbsLoc')

	exonLen         =  os.path.join(workfloder,'exonLen')	 
	exonjunfile     =  os.path.join(workfloder,'GeneExonFile')	 
	spandExonsLen   =  os.path.join(workfloder,'exonLen')

	
	static_exonJunc.staticGeneNewExon(GeneExonSplit,geneMapout,geneExonout)
	getIsoLength.getIsoLen(GeneFile,GeneExonSplit,targetDir)#get isoform length
	geneExonLen.staticGeneNewExon(GeneExonSplit,exonjunfile,spandExonsLen,ReadLen)	

	plusGeneName.plusGeneName(GeneFile,readFileList,targetDir,readType)#covert bowtie
	
	DictIsofGene = {}
	gene_list=[]
	isoNo_list=[]
	length_list=[]
	isoName_list=[]
	f = open(GeneFile,'r')
	for line in f:
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]
		IsoNo = int(line[1])
		IsoNameList=  line[3:]
		for i in xrange(IsoNo):
			DictIsofGene[IsoNameList[i]]=gene

		
		gene_list.append(line[0])
		isoNo_list.append(int(line[1]))
		length_list.append(int(line[2]))
		isoName_list.append(line[3:])
	f.close()


    ###########CalculateProbability and ExtractGeneData  
	for i in range(len(Sample)):                           
		for j in range(len(rep)):
			InputPath = os.path.join(workfloder,'readInput')
			InputFile = os.path.join(InputPath,'Lane'+str(i)+'.plusGene')
			FileLabel = os.path.join(workfloder,'NLDMseq.Gene.Label.Isoform')
			OutputFile= os.path.join(InputPath,'Lane'+str(i))

			seq_depth=Extract_BowtieData.CalculateProbability(InputFile,OutputFile,readType)

			InputFile = OutputFile

			OutputPath= os.path.join(GeneData_Path,lanFile[j])		

			Extract_BowtieData.ExtractGeneData(InputFile,OutputPath,GeneFile,FileLabel,readType)

                          	   
	##################CalculateAbsoluteLocation	
	for l in range(len(rep)):
		CalLocFileInPath  = os.path.join(GeneData_Path,lanFile[l])
		CalLocFileOutPath = os.path.join(GeneAbsLoc_path,lanFile[l])
		Extract_CalcAbsLoc.CalculateAbsoluteLocation(GeneFile,GeneExonSplit,CalLocFileInPath,CalLocFileOutPath,readType,ReadLen)



    	  
	targetDir=AnnotationDir
	LocInputPath   = os.path.join(targetDir,'workfloder','ModelMultiGene_AbsLoc')
	exonInputFile  = os.path.join(targetDir,'workfloder','GeneExonFile')
	exonLenFile    = os.path.join(targetDir,'workfloder','exonLen')
	DataOutputPath = os.path.join(targetDir,'workfloder','ModelMultiGene_Data')
	NormDataOutput = os.path.join(targetDir,'workfloder','ModelMultiGene_NormData')	 

	#function pp_work parament list and short name:
	# 	sub_gene_list	----->sgl
	#  sub_length_list	----->sll
	#   sub_isoNo_list	----->siNol
	# sub_isoName_list	----->siNamel
	#     LocInputFile	----->LocIF
	#    exonInputFile	----->exonInF
	#      exonLenFile	----->exonLenF
	#   DataOutputPath	----->DOPath
	#   NormDataOutput	----->NDOut
	def pp_work(sgl,sll,siNol,siNamel,LocIF,exonInF,exonLenF,DOPath,NDOut):
		for i in range(len(sgl)):
			gene=sgl[i];
			length=sll[i];
			isoNo=siNol[i];
			isoName=siNamel[i];
			static_readOnExon.ModelMultiGeneDataScale(gene,length,isoNo,isoName,LocIF,exonInF, \
															exonLenF,DOPath,NDOut);
	


	for i in range(len(Sample)):                           
		    for j in range(len(rep)):
				LocInputFile = os.path.join(LocInputPath,lanFile[j])
				block_size=5000
				start=0
				end=len(gene_list)
				job_server = pp.Server()
				job_server.get_ncpus()
				jobs=[]
				endi=0;
				starti=0;
				while endi<end:
					endt=starti+block_size;
					endi=min(endt,end);
					jobs.append(job_server.submit(pp_work,(gene_list[starti:endi], \
															length_list[starti:endi],\
															isoNo_list[starti:endi],\
															isoName_list[starti:endi],\
															LocInputFile,\
															exonInputFile,\
															exonLenFile,\
															DataOutputPath,\
															NormDataOutput),\
															globals=globals()))
					starti=starti+block_size;
				for job in jobs:
					job()
				job_server.destroy()
	return (seq_depth)