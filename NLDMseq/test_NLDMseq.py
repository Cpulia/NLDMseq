#!/usr/bin/env python
#-*- coding:utf-8 -*-

import os 
import os.path
import shutil
import ctypes

if not os.path.exists("libgene_expre.so"):
        print "Build software filed.Please rebuild again use the install.sh script"
if os.path.exists("./workfloder"):
	shutil.rmtree("./workfloder")

os.mkdir("./workfloder")
if not os.path.exists("./workfloder"):
        print "Create workfloder failed .Please get the access control"
else:
        so=ctypes.CDLL("./libgene_expre.so")
        so.gene_expre("./TEST_DATA","test_gene",100,"./workfloder")
        if os.path.exists("./workfloder/expression_log"):
                fp=open("./workfloder/expression_log")
                expre_data=fp.readline()
                gene_name=expre_data.split()
                if not (gene_name[0]=="test_gene"):
                        print "Test filed .Please try again."
                shutil.rmtree("./workfloder")
		print "Test passed! The software has been installed correctly"
