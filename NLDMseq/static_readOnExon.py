#!/usr/bin/env python
# -*- coding: utf-8 -*-


from __future__ import division
import os
import os.path
 
 




def ModelMultiGeneDataScale(GeneName,GeneLen,IsoNo,IsoName,LocInputFile,ExonlistFile,ExonLenPath,DataOutput,NormDataout):
        
        ReadCout_onExon={}
        f_exon=open(ExonlistFile+'/'+GeneName,'r');
        exonlist=f_exon.readline();
        exonlist=exonlist.rstrip()
        exonlist=exonlist.split('\t')
        f_exon.close()
        f_exon=open(ExonLenPath+'/'+GeneName,'r')
        exonLen=f_exon.readline()
        exonLen=exonLen.rstrip()
        exonLen=exonLen.split('\t')
        f_exon.close()
        
        readCount = [0 for i in range(len(exonlist))]
        K=0
          
        targetFile = os.path.join(LocInputFile, GeneName);
        if os.path.isfile(targetFile):
          
           f_in = open(LocInputFile+'/'+ GeneName,'r');
          # print GeneName

           dic_read = dict()
           for line in f_in:
               line = line.rstrip();
               line = line.split('\t');
               readid = line[0];

               if readid not in dic_read:
                  dic_read[readid] = 1;
               else:
                  dic_read[readid] = dic_read[readid] + 1;

           f_in.close();
           
           LossRead = 0 
           f_in = open(LocInputFile+'/'+ GeneName,'r');
           while True:
               line = f_in.readline();
               if not line: break;
               line = line.rstrip();
             
               line = line.split('\t');
            # if int(len(line))==4:
               readid = line[0];
              # print line[1]
               iso_nm = line[1];
               #print iso_nm
               abs_loc = int(line[2]);
               exonIndex = line[3];
               prob = float(line[4])
               abs_ary = [0 for i in range(int(IsoNo))];  
               iso_flag = [0 for i in range(int(IsoNo))];
               exon_ary = ['1' for i in range(int(IsoNo))]
               ind = IsoName.index(iso_nm);
               abs_ary[ind] = abs_loc
               iso_flag[ind] = 1;
               exon_ary[ind]= exonIndex;

               if readid in dic_read:
                  count = dic_read[readid];
               flag = 1;
               while flag != count:
                     flag = flag + 1;
                     info = f_in.readline();
                     info = info.rstrip();
                     info = info.split('\t');
                     if len(info)==5:
              
		                 info_iso = info[1];
		                 info_loc = int(info[2]);
		                 inf_exon= info[3] 
		
		                 ind = IsoName.index(info_iso);
		                 abs_ary[ind] = info_loc;
		                 iso_flag[ind] = 1;
		                 exon_ary[ind] = inf_exon;
		              
               ind = iso_flag.index(1);
               ind_loc = abs_ary[ind];
               ind_exon = exon_ary[ind]
              
               loc_temp = []
               exon_temp = []
               for i in range(IsoNo):
                       if iso_flag[i] == 1:
                          loc_exon=str(exon_ary[i])+':'+str(abs_ary[i])
                          loc_temp.append(loc_exon)
                          
                       else:
                          pass
#                   print loc_temp, len(loc_temp)
               filer_Dic={}
               for i in range(len(loc_temp)):
                       if loc_temp[i] not in filer_Dic:
                          filer_Dic[loc_temp[i]]=1;
                       else:
                          filer_Dic[loc_temp[i]]=filer_Dic[loc_temp[i]]+1
                       

           #print str(len(filer_Dic))+'\t'+str(len(loc_temp))
               
               #print filer_Dic
               #print loc_temp
               read_exons=filer_Dic.keys();
                   
               
               for i in range(len(filer_Dic)):#
                   exon=read_exons[i];#
                   exon=exon.split(':');
                   exon=exon[0]
                   
                   if exon in exonlist: 
                     readCount[exonlist.index(exon)]=readCount[exonlist.index(exon)]+1#GGGGG
                     #readCount[exonlist.index(exon)]=readCount[exonlist.index(exon)]+prob/len(filer_Dic);
                   else:
                     #readCount[0]=readCount[0]+prob/len(filer_Dic)
                     readCount[0]=readCount[0]+1#GGGGG

           #print str(readCount)
           f_out=open(DataOutput+'/'+GeneName,'a');
           for L in range(len(readCount)):
            #f_out.write(str(L+1)+':'+str(GeneBlockDat[K][L])+'\t')
               f_out.write(str(L+1)+':'+str(int(readCount[L]))+'\t')
           f_out.write('\n') 
           #print "readcout:"+str(readCount)
          # print  "exonLen:"+str(exonLen)
           f_out=open(NormDataout+'/'+GeneName,'a')
           
         
           for L in range(len(readCount)):
                a=int(exonLen[L]);
                if a>0:
                   
                   #f_out.write(str(L+1)+':'+str(  (int(readCount[L])/a+0.001)*1000 )+'\t')
                   f_out.write(str(L+1)+':'+str(  ((readCount[L])/a+0.001))+'\t')#GGGGGG
                
                else:
                   f_out.write(str(L+1)+':'+str(1)+'\t')
           f_out.write('\n')
ReadLen=35
