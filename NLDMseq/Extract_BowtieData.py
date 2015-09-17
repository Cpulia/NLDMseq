#!/usr/bin/env  python
# -*- coding: utf-8 -*-
from __future__ import division  


import pp
import os, os.path
import random as rd
import numpy as np
import gc



def CalculateProbability(InputFile,OutputFile,ReadType):

    DictRead = {}

    f = open(InputFile,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        read = line[0]
        if ReadType == 'Single':
            gene = line[3]
        if ReadType == 'Paired':
            gene = line[4] 

        if read not in DictRead:
            genelist = [gene]
            DictRead[read] = genelist
        else:
            # genelist = DictRead[read]
            if gene in genelist:
                pass
            else:
                genelist.append(gene)
                DictRead[read] = genelist
    f.close()

    GeneUniqRead = {}

    for key in DictRead.iterkeys():
        read = key
        genelist = DictRead[key]

        if len(genelist) == 1:
            gene = genelist[0]
            if gene not in GeneUniqRead:
               GeneUniqRead[gene] = 1
            else:
               GeneUniqRead[gene] = GeneUniqRead[gene] + 1
        else:
            pass


    for key in DictRead.iterkeys():
        read = key
        genelist = DictRead[key]
        geneno = len(genelist)

        if geneno == 1:
            geneRate = [1.0]
            DictRead[read] = [genelist,geneRate]
        else:
            countlist = []
            for i in xrange(geneno):
                gene = genelist[i]
                if gene in GeneUniqRead:
                    count = GeneUniqRead[gene]
                else:
                    count = 0 
                countlist.append(count)

            if sum(countlist) ==0:                   # 
                rd.seed(123456)
                rdno = rd.randint(0,geneno-1)
                geneRate = [0.0 for j in xrange(geneno)]
                geneRate[rdno] = 1.0
            else:
                MaxGene = max(countlist)
                MaxInd = countlist.index(MaxGene)
                geneRate = [0.0 for j in xrange(geneno)]
                geneRate[MaxInd] = 1.0

            DictRead[read] = [genelist,geneRate]          

    

    f = open(InputFile,'r')
    fo = open(OutputFile+'.Prob','w')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        
        if ReadType == 'Single':
            read, loc, trans, gene = line
        if ReadType == 'Paired':
            read, loc1, loc2, trans, gene = line



        genelist, geneRate = DictRead[read]
        geneInd = genelist.index(gene)
        readRate = geneRate[geneInd]
        
        if ReadType == 'Single':
            fo.write(read+'\t'+gene+'\t'+trans+'\t'+loc+'\t'+str(readRate)+'\n')
        if ReadType == 'Paired':
            fo.write(read+'\t'+gene+'\t'+trans+'\t'+loc1+'\t'+str(loc2)+'\t'+str(readRate)+'\n')       
    
    f.close()
    fo.close()

    flab = open(OutputFile+'.Read','w')
    readNo = 0
    for read in DictRead.iterkeys():
        flab.write(read+'\t'+str(readNo)+'\n')
        readNo = readNo+1
    flab.close()


    del DictRead
    gc.collect()

    del GeneUniqRead
    gc.collect()


    return readNo


def ExtractGeneData(InputFile,OutputPath,ObjectGene,FileLabel,ReadType):


    DictTransLab = {}
    DictGeneLab = {}
    f = open(FileLabel,'r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        gene,lab, trans = line

        if gene not in DictGeneLab:
            DictGeneLab[gene] = [trans]
        else:
            temp  = DictGeneLab[gene]
            temp.append(trans)
            DictGeneLab[gene] =temp        
        DictTransLab[trans] = int(lab)
    f.close()

    DictReadLab = {}
    DictReadLabC = {}
    f = open(InputFile+'.Read','r')
    for line in f:
        line = line.rstrip()
        line = line.split('\t')
        read, lab = line     
        DictReadLab[read] = lab
        DictReadLabC[lab] = read
    f.close()




    DictGeneData = {}    
    if ReadType == 'Single':
        lineNo = 1
        flag = 0
        fi = open(InputFile+'.Prob','r') 
        for line in fi:
            line = line.rstrip()
            line = line.split('\t')
            read,gene,trans,loc,prob = line

            readlab = DictReadLab[read]
            translab = DictTransLab[trans]


            if gene not in DictGeneData:
                readlist = [readlab]
                problist = [prob]
                translist = [translab]
                loclist = [loc]
                DictGeneData[gene]= readlist,problist,translist,loclist
            else:
                readlist,problist,translist,loclist = DictGeneData[gene]
                readlist.append(readlab)
                problist.append(prob)
                loclist.append(loc)
                translist.append(translab)                
                DictGeneData[gene]= readlist,problist,translist,loclist

            if lineNo%10000000==0:
                flag = 1
            else:
                flag = 0

            lineNo = lineNo+1

            if flag ==1:
                print '.',
                f = open(ObjectGene,'r')
                for line in f:
                    line = line.rstrip()
                    line = line.split('\t')
                    gene = line[0].rstrip()           

                    if gene in DictGeneData:
                        readlist,problist,translist,loclist = DictGeneData[gene]
                        transName = DictGeneLab[gene]

                        geneFile = OutputPath+'/'+gene
                        if os.path.isfile(geneFile):
                            fo = open(geneFile,'a')
                        else:
                            fo = open(geneFile,'w')

                        for i in range(len(readlist)):
                            read = DictReadLabC[readlist[i]]
                            trans = transName[translist[i]-1]
                            fo.write(read+'\t'+trans+'\t'+loclist[i]+'\t'+problist[i]+'\n')
                        fo.close()
                    else:
                        pass                    
                f.close()  
                DictGeneData = {}
                gc.collect()
                flag = 0
            else:  
                pass     
        fi.close()


        f = open(ObjectGene,'r')
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            gene = line[0].rstrip()           

            if gene in DictGeneData:
                readlist,problist,translist,loclist = DictGeneData[gene]
                transName = DictGeneLab[gene]

                geneFile = OutputPath+'/'+gene
                if os.path.isfile(geneFile):
                    fo = open(geneFile,'a')
                else:
                    fo = open(geneFile,'w')

                for i in range(len(readlist)):
                    read = DictReadLabC[readlist[i]]
                    trans = transName[translist[i]-1]
                    fo.write(read+'\t'+trans+'\t'+loclist[i]+'\t'+problist[i]+'\n')
                fo.close()
            else:
                pass                    
        f.close()  

        print lineNo
    elif ReadType == 'Paired':
        lineNo = 1
        flag = 0       
        fi = open(InputFile+'.Prob','r')         
        for line in fi:
            line = line.rstrip()
            line = line.split('\t')
            read,gene,trans,loc1,loc2,prob = line
    
            readlab = DictReadLab[read]
            translab = DictTransLab[trans]
            
            if gene not in DictGeneData:
                readlist = [readlab]
                problist = [prob]
                translist = [translab]
                loc1list = [loc1]
                loc2list = [loc2]
                DictGeneData[gene]= readlist,problist,translist,loc1list,loc2list
            else:
                readlist,problist,translist,loc1list,loc2list = DictGeneData[gene]
                readlist.append(readlab)
                problist.append(prob)
                loc1list.append(loc1)
                loc2list.append(loc2)
                translist.append(translab)                
                DictGeneData[gene]= readlist,problist,translist,loc1list,loc2list

            if lineNo%10000000==0:
                flag = 1
            else:
                flag = 0

            lineNo = lineNo+1
            
            if flag ==1:
                f = open(ObjectGene,'r')
                for line in f:
                    line = line.rstrip()
                    line = line.split('\t')
                    gene = line[0].rstrip()           

                    if gene in DictGeneData:
                        readlist,problist,translist,loc1list,loc2list = DictGeneData[gene]
                        transName = DictGeneLab[gene]

                        geneFile = OutputPath+'/'+gene
                        if os.path.isfile(geneFile):
                            fo = open(geneFile,'a')
                        else:
                            fo = open(geneFile,'w')

                        for i in range(len(readlist)):
                            read = DictReadLabC[readlist[i]]
                            trans = transName[translist[i]-1]
                            fo.write(read+'\t'+trans+'\t'+loc1list[i]+'\t'+loc2list[i]+'\t'+problist[i]+'\n')
                        fo.close()
                    else:
                        pass                    
                f.close()  
                DictGeneData = {}
                gc.collect()
                flag = 0
            else:  
                pass 
        fi.close()


        f = open(ObjectGene,'r')
        for line in f:
            line = line.rstrip()
            line = line.split('\t')
            gene = line[0].rstrip()           

            if gene in DictGeneData:
                readlist,problist,translist,loc1list,loc2list = DictGeneData[gene]
                transName = DictGeneLab[gene]

                geneFile = OutputPath+'/'+gene
                if os.path.isfile(geneFile):
                    fo = open(geneFile,'a')
                else:
                    fo = open(geneFile,'w')

                for i in range(len(readlist)):
                    read = DictReadLabC[readlist[i]]
                    trans = transName[translist[i]-1]
                    fo.write(read+'\t'+trans+'\t'+loc1list[i]+'\t'+loc2list[i]+'\t'+problist[i]+'\n')
                fo.close()
            else:
                pass                    
        f.close()                     


    else:
        pass


    del DictTransLab
    del DictGeneLab
    gc.collect()
    del DictReadLab
    gc.collect()
    del DictReadLabC
    gc.collect()   
    del DictGeneData
    gc.collect()

