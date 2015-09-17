#!/usr/bin/env python
#coding=utf-8


import os
import os.path
import sys
import getopt
import shutil
import time

 
#my own lib
#Annoation
import analyseAnnotation
# # pre-processing
import Extract_GeneData
#LDA 
import lda

DEBUG=True  
 
if DEBUG == True: 
        import pdb 
########################################
#function:give a sample help list when the software does not work 
########################################
def usage():
    print sys.argv[0] +'   -a/--AnnotationFile: The annotation includes the gene and isoform information. eg: refGene, knownGene, and ensGene, which all can be download UCSC website. If you use the Ensembl dataset, you may use the .gtf file.'
    print sys.argv[0] +'   -t/--AnnotationType   <int>    supported four annotation types: refGene: 1, ensGene: 2, knownGene: 3 and Ensembl: 4'
    print sys.argv[0] +'   -i/--Input: Alignment File(s),align_reads.output(s) in a condition,use comma separating the input files if more than one files are provided.(e.g,align_reads1.output,align_reads2.ouput)'
    print sys.argv[0] +'   -o/--OutputPath: All output files are in this direction  '
    print sys.argv[0] +'   -s/--SingleEnd: Input data are Single-end reads alignment. (Default: off)'
    print sys.argv[0] +'   -h/--Help #get help info'
#function:check the refference files ,input files and output path of the software
def check_path(ref_file,input_path,output_path):
    if( not os.path.exists(ref_file)):
        print sys.argv[0]+'    Error:'+'The annotation file does not exist.'
        usage()
        sys.exit()
    for item in input_path:
        if( not os.path.exists(item)):
            print sys.argv[0]+'    Error:'+'The input files does not exist.'
            usage()
            sys.exit()

#function:get the readlen from read data_exist
def get_readlen(read_path,read_format):
    readlen=0
    if(read_format=='single_end'):
        fp=open(read_path[0],'r')
        line=fp.readline()
    #                line=line.rsplit()
        line=line.split('\t')
        readlen=len(line[10])  
        print("the read length of input file is %s\n" %str(readlen))
    else:
        fp=open(read_path[0],'r')
        line=fp.readline()
    #                line=line.rsplit()
        line=line.split('\t')
        readlen=len(line[10])  
        print("one of read tail length of input file is %s\n" %str(readlen))         
    return(readlen)                 

########################################
#function:create the work path for all the 
#input:work_space,readFileList
#output:readlen
 ########################################
def CreateRunFileName(work_space,readFileList):
    laneFile=[]
    for i in range(len(readFileList)):
        laneFile.append('Lane'+str(i+1))
    
    workfloder      = os.path.join(work_space, 'workfloder')
    GeneData_Path   = os.path.join(workfloder,'ExtractMultiGeneData')
    GeneAbsLoc_path = os.path.join(workfloder,'ModelMultiGene_AbsLoc')

    os.mkdir(GeneData_Path)
    os.mkdir(GeneAbsLoc_path)

    os.mkdir(os.path.join(workfloder,'ModelMultiGene_Data'))
    os.mkdir(os.path.join(workfloder,'ModelMultiGene_NormData'))

    for lan in laneFile:
        sub_GeneData_path   = os.path.join(GeneData_Path,lan);
        sub_GeneAbsLoc_path = os.path.join(GeneAbsLoc_path,lan)

        os.mkdir(sub_GeneData_path);
        os.mkdir(sub_GeneAbsLoc_path)

    readInput_path = os.path.join(workfloder,'readInput')
    mapFile        = os.path.join(workfloder,'ModelMultiGene_Map')
    GeneExonFile   = os.path.join(workfloder,'GeneExonFile')     
    isoLenFile     = os.path.join(workfloder,'isoLen')
    exonLen        = os.path.join(workfloder,'exonLen')

    os.mkdir(readInput_path)
    os.mkdir(mapFile)
    os.mkdir(GeneExonFile)
    os.mkdir(isoLenFile)
    os.mkdir(exonLen)



ref_exist=False
ref_type_exist=False
input_exist=False
out_exist=False
readType='Paired'

# get the parameter from the cmd line 
opts, args = getopt.getopt(sys.argv[1:], "hsa:i:o:t:", ["Help","PairedEnd", "AnnotationFile","Input", "OutputName","AnnotationType"])
for op, value in opts:
    if op == '-a' or op == '--AnnotationFile':
        ref_file=value
        ref_exist=True
    elif op in ["-t","--AnnotationType"]:
        typeflag = value
        ref_type_exist=True
    elif op == '-i' or op == '--Input':
        input_path = value
        input_path=input_path.split(',')
        input_exist=True
    elif op == '-o' or op == '--OutputName':
        output_path = value
        out_exist=True
    elif op == '-s' or op == '--SingleEnd':
        readType='Single'
    elif op == '-h'or op == '--Help':
        usage()
        sys.exit()
# check the cmd                
data_exist=ref_exist and input_exist and out_exist and ref_type_exist
if(not data_exist):
   usage()
   sys.exit()
       
start_time=time.strftime('%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
  
if DEBUG== True:
    pdb.set_trace(); 

check_path(ref_file,input_path,output_path)
output_path=os.path.abspath(output_path)
work_path= output_path
work_floder_path=os.path.join(work_path,'workfloder')


if( os.path.exists(work_floder_path)):
    shutil.rmtree(work_floder_path)      
os.makedirs(work_floder_path)
CreateRunFileName(work_path,input_path)


read_length=get_readlen(input_path,readType)

AnnoType = {'1':'refGene','2':'ensGene','3':'knownGene','4':'Ensembl'}
annotation_type = AnnoType[typeflag]
annotation_file=ref_file
prefix_name='NLDMseq'
outputDir=work_floder_path

analyseAnnotation.CovertAnnotation(annotation_type,annotation_file,prefix_name,outputDir)
analyseAnnotation.AnalyseAnnotation(annotation_file,prefix_name,outputDir)      



seq_depth=Extract_GeneData.startProcess(input_path,work_path,read_length,readType)
print seq_depth;
print 'Calculate... ' 
lda.lda_expre(work_floder_path,seq_depth,output_path)
# lda.removeFileInFirstDir(work_floder_path)
# shutil.rmtree(work_floder_path)    
stop_time=time.strftime('%Y-%m-%d-%H-%M-%S',time.localtime(time.time()))
print('start time:%s\n' %start_time)
print('stop time:%s\n' %stop_time)  
