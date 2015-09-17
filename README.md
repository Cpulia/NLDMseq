#NLDMseq
NLDMseq is the software for expression calculation at both gene and isoform levels from RNA-seq data. The program calculates expression values using the alignment from Bowtie 2.
#What is NLDMseq?

NLDMseq is the software for expression calculation at both gene and isoform levels from RNA-seq data given a reference transcriptome. The program calculates expression values using the alignment from Bowtie 2.

The software is an open-source software. You can redistribute it and/or modify it under the terms of the GNU General Public License (GPL) as published by the Free Software Foundation.

#How NLDMseq works?

##Features

>* NLDMseq adopts latent variables to model the unknown isoforms, from which reads originate, and the underlying percentage of multiple spliced variants.

>* NLDMseq provides an approach to accurately estimate gene and isoform expression from RNA-Seq data by modeling the isoform- and exon-specific read sequencing biases.
>* NLDMseq provides the measurement uncertainty associated with the estimated gene and isoform expression.

#Contact
Users can ask technical questions by sending emails to Xuejun Liu (xuejun.liu@nuaa.edu.cn).

#Installation

##Installation requirements:

>* Operating system :
	Linux( Ubuntu13.10 or Fedora 17…)
	Mac OS X (Mac OS X 10.8 5 or higher)
>* Python 2.7 
>* GUN Scientific Library GCC 4.8.1 or higher(http://www.gnu.org/software/gsl/)
>* Python package: parallel python (pp)(http://www.parallelpython.com/) 
>* Bowtie2

###NOTE: 
 PP is a python module which provides mechanism for parallel execution of python code on SMP(systems with multiple processors or cores) and clusters (computers connected via network).It is light, easy to install and integrate with other python software.PP is an open source and cross-platform module written in pure python .

The software depends on the free and open-source software, the GNU Scientific Library (GSL) (http://www.gnu.org/software/gsl/), so the GSL needs to be installed on the user’s system. NLDMseq requires the GSL installed in /usr/local, which is the default location. The GSL can be compiled by the user. Users who are compiling NLDMseq from the source code should install GSL in the standard location (/usr/local). If you are not sure if GSL is already installed, at the Terminal prompt $ type:`$gsl-config --cflags --libs-without-cblas`
If GSL is installed, the command above should return the following information:
```shell
-I/sw/include
-L/sw/lib -lgsl –lm
```
GSL can be found in the gsl subdirectory on your nearest GNU mirror ( http://ftpmirror.gnu.org/gsl/) or the Main GNU ftp site (ftp://ftp.gnu.org/gnu/gsl/). The users can dowmload gsl-1.6.tar.gz or a higher version for NLDMseq. Please follow the instructions in the included file “INSTALL” to guide  the installation of this library.
We recommend using the Linux operating system. The test shows that it will cost much more time on a Mac OS X machine.

##Installation

>* First, go to the location of NLDMseq

```
    $cd ~
```

>* Second, install the software by the install.sh shell script.

```
$./INSTALL.sh
```

NOTE:The user may have to give the script execute permissions by chmod u+x INSTALL.sh.This script will use the gcc complier to build C codes in the software.
>* Test the software
```shell
$python test_NLDMseq
```

To test this software, the user can run the python script test_NLDMseq.py.This will check whether the software has been built and installed correctly. It will calculate the expression of genes and isoforms using the data from the TEST_DATA folder.

##Running NLDMseq
>* Step 1. Aligning sequenced reads with Bowtie 2
The following commands are used to align sequenced reads to a reference transcriptome.
```
$ bowtie2-build -f ref_transcript. Fasta ref_transcript.index
$ bowtie2 –t –f -k 20 -p 4--no-sq --no-head -x ref_transcript.index raw_data.fasta -S align_reads.output
```
If the paired-end reads are processed, the Bowtie command should be listed as
```
$ bowtie2 –t–f -k 20 -p 4 --no-sq --no-head -x ensGene.ref_transcript.index -1 raw_data.fasta -2 raw_data.fasta -S align_reads.output
```

The above transcriptome reference sequence can be downloaded from UCSC or Ensembl website.

>* Step 2. Calculate expression values with NLDMseq

The software uses the alignment SAM file from Bowtie 2 as input data.

```
$ python NLDMseq-calulation.py –a AnnotationFile.txt  -t AnnotationFile Type –i align_file(s) –o  expression_path
```

Options:

>* -a/--AnnotationFile: The annotation includes gene and isoform information. (e.g: refGene, knownGene, and ensGene, which can all be downloaded from UCSC website.)
>* -t/--AnnotationType <int>: an integer to select the type of annotation. refGene:1, ensGen:2, knownGene:3 and Ensmbl: 4. For detailed interpretation, please click here.
>* -i/--Input: Alignment file(s), which is Bowtie 2 output.More than one output files for multiple lanes under a condition.use comma for the separation(e.g, align_reads1.output, align_reads2.ouput).
>* -o/--OutputPath All output files and temp files in the process of calculation are produced in this given path. (E.g. two output files will be produced, named expression_gene.txt and expression_isofrom.txt)
>* -s/--Single-end: Optional. Input data are alignment for single-end reads. (Default: off)
>* -h/-Help: Show help information.

Output files:
The NLDMseq produces gene and isoform expression output files, which are put in the location by the path given in the option -o/--OutputPaht


Table 1: Description of output files:

|Column Nr |Name |Description|
|:--------:|:----------:|:--------:|
|1|gene/isoform Name|expression_gene.txt:Name of a gen                 expression_isoform.txt:Name of gene and isoform|
|2| FPKM| FPKM expression value of gene/isoform.|
3 |VarianceFPKM|Variance value of gene/isoform.|

#Example
Here, we use a simple example to show the usage of NLDMseq. The alignment files from Bowtie 2 and the annotation file are also supplied.

>* Annotation file: ensGene.txt
>* Alignment from Bowtie2: SE_read_example.sam (single-end)
>* Alignment from Bowtie2:PE_read_example.sam (paired-end)

Since the Bowtie2 output has been supplied, so you can skip step1 and just run the following command.

For single-end data:
```
$ python NLDMseq-calulation.py –a ../ensGene.txt  -t 2 -i ../SE_read_example –o ../se_out -s 
```
After running the above command, you will obtain two output files of this single-end data, expression_gene and expression_isoform under se_out path.

For paired-end data :
```
$ python NLDMseq-calulation.py –a ../ ensGene.txt –t 2 –i ../PE_read_example –o../pe_out
```
After running the above command, you will obtain two output files of this paired-end data, expression_gene and expression_isoform under se_out path.

#Authors
Xuejun Liu (xuejun.liu@nuaa.edu.cn).College of Computer Science and Technology, Nanjing University of Aeronautics and Astronautics, 29Jiangjun Rd., Jiangning District, 211106Nanjing China.
We use the C source code from Daichi Mochihashi (ATR Spoken Language Communication Research Laboratories, Kyoto, Japan daichi.mochihashi@atr.jp)
#Related software

Bowtie2: for alignment of RNA-Seq reads to transcriptome (optionally with the precomputed set of splice junctions).

