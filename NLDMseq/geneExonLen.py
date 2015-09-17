#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import os
import os.path
import string

def staticGeneNewExon(InputFile,exon_juncNum,output,ReadLen):
    print 'Begin statistic exon lengh...'
    f_in = open(InputFile)

    for line in f_in:
      line = line.rstrip()
      line = line.split()
      gene = line[0]
     
      iso_no = line[1]
      if  iso_no.isdigit():
 
          exon_no=int(line[2])
          exonLenList=[]
          lenTemp=line[4].split(',')
          exonLenList.append(int(lenTemp[0]))
          for i in range(exon_no-1):
            exonLenList.append(int(lenTemp[i+1])-int(lenTemp[i]))
           
           
         # print exonLenList
          if os.path.isfile(os.path.join(exon_juncNum,gene)): 
              f_gene = open(exon_juncNum+'/'+gene,'r');
              line=f_gene.readline();
              line=line.rstrip()
              line=line.split('\t')
             # print line
              f_gene.close()

              for i in range(len(line)-exon_no):
                 exonLenList.append(ReadLen-1)
          #print exonLenList
          f_open=open(output+'/'+gene,'w')
          for i in range(len(exonLenList)):
              f_open.write(str(exonLenList[i])+'\t')
          f_open.close()
    f_in.close() 

 
