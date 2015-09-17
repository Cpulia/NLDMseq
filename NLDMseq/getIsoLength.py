#! /usr/bin/env python

import os
import os.path
def getIsoLen(refFlat_Check_GeneUnionExonLenFile,refFlat_Check_MultiGeneExonSplitFile,targetDir):
	print 'Begin get isoform length...'

	targetDir  = os.path.join(targetDir,'workfloder')
	isoLenFile = os.path.join(targetDir,'isoLen')
	if not os.path.exists(isoLenFile):
		 os.mkdir(isoLenFile)


	f_in = open(refFlat_Check_MultiGeneExonSplitFile,'r');

	DictIsofLength = {}

	for line in f_in:
		line = line.rstrip();
		line = line.split('\t');
		isofName = line[1]
		exonNo = int(line[2])

		LenTemp= line[4].split(',')
		LenTemp.pop()
		IsoLen = LenTemp[exonNo-1]

		if isofName not in DictIsofLength:
			DictIsofLength[isofName] = IsoLen
		else:
			pass

	f_in.close();




	f_isoin = open(refFlat_Check_GeneUnionExonLenFile,'r')
	for line in f_isoin:
		line = line.rstrip();
		line = line.split('\t')

		geneName = line[0]
		isoNo = int(line[1])

		isoNameList = line[3:]

		f_out = open(isoLenFile+'/'+line[0],'w')
		for i in xrange(isoNo):
			f_out.write(DictIsofLength[isoNameList[i]]+'\n')
		f_out.close()


	f_isoin.close()

	del DictIsofLength
