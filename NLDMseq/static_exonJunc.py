#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import os
import os.path
def staticGeneNewExon(InputFile,MapOutputFile,geneExonFile):
    
    
    
    print 'Begin creating preudo-exon for every gene'
    f_in = open(InputFile)
    line = f_in.readline();
    i = 0 ;
    while(line):
		 
		
		line = line.rstrip()
		line = line.split('\t')
		gene = line[0]
		iso_no = int(line[1])
		flag=0
		isoExonList = [[]for n in range(iso_no)]
		isoBetaMatrix = [[]for n in range(iso_no)]
		if (flag==0):               
			geneExonList = line[3].split(',')
			flag=1;
			for j in range(iso_no):
				line = f_in.readline()
				line = line.split('\t')
				isoExonList[j] = line[3].split(',');
				 
				for k in range(len(isoExonList[j])-1):
					conjunc = str(isoExonList[j][k])+'-'+str(isoExonList[j][k+1])
					isoExonList[j].append(conjunc);#
					if conjunc not in geneExonList:
						geneExonList.append(conjunc);#
					else:
						pass;
				 
			#print gene;
			#print geneExonList;   
			
			for j in range(len(geneExonList)):#
				for m in range(len(isoExonList)):
					if geneExonList[j] not in isoExonList[m]:
						isoBetaMatrix[m].append(0)
					else:
						isoBetaMatrix[m].append(1)
			if (iso_no>=1):
			 
                         f_out = open(MapOutputFile+'/'+ gene,'w');
                         f_exon = open(geneExonFile+'/'+gene,'w')
                         for k in range(len(geneExonList)):
                               f_exon.write(str(geneExonList[k])+'\t')
				
                         for K in range(iso_no):
					for L in range(len(geneExonList)):
						f_out.write(str(isoBetaMatrix[K][L])+'\t')
					f_out.write('\n')
				
        
			
			
			
			i = i+1#
		 
			line=f_in.readline();
    f_out.close()			
    f_exon.close()
    print i;	 
    f_in.close();
	
 
 




            
