/*
    lda.c
    Latent Dirichlet Allocation, main driver.
    (c) 2004 Daichi Mochihashi, All Rights Reserved.
    $Id: lda.c,v 1.11 2004/11/10 04:23:06 dmochiha Exp $

*/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include<string.h>
#include<math.h>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "lda.h"
#include "learn.h"
#include "writer.h"
#include "feature.h"
#include "dmatrix.h"
#include "util.h"
#include "read_beta.h"
#include  "gene_expression.h"

#define SAMPLE_CNT 1000
/*int main()
{
    gene_expre("/home/cpulia/Desktop/NLDM_test/SE/workfloder","ENSG00000136868",100,"/home/cpulia/gene_exp");
    return(0);
}*/
int
gene_expre (char* data_path,char* gene_name,int SEQ_DEPTH,char *output_path)
{
	document *data=NULL;
	document *sum_exon=NULL;//creat the space for the sum of read tunnel
	double *alpha;
	double *short_alpha;
	double *theta;
	double **beta;
    double ** isoform;
    double *iso_expre;
    double *iso_var;
    double *isolen;
	int nlex, dlenmax,lda_nclass;
	int nclass     = CLASS_DEFAULT;		// default in lda.h

    int beta_size[2];//line and row of beta
    /*input data path*/
    char *gene_data_path=(char*)calloc(256,sizeof(char));
    gene_data_path=strcat( gene_data_path,data_path);
    gene_data_path=strcat( gene_data_path,"/ModelMultiGene_Data/");
    gene_data_path=strcat( gene_data_path,gene_name);

    char *gene_norm_data_path=(char*)calloc(256,sizeof(char));
    gene_norm_data_path=strcat(gene_norm_data_path,data_path);
    gene_norm_data_path=strcat(gene_norm_data_path,"/ModelMultiGene_NormData/");
    gene_norm_data_path=strcat(gene_norm_data_path,gene_name);

    char *isoLen_path=(char*)calloc(256,sizeof(char));
    isoLen_path=strcat(isoLen_path,data_path);
    isoLen_path=strcat(isoLen_path,"/isoLen/");
    isoLen_path=strcat(isoLen_path,gene_name);

    char *map_path=(char*)calloc(256,sizeof(char));
    map_path=strcat(map_path,data_path);
    map_path=strcat(map_path,"/ModelMultiGene_Map/");
    map_path=strcat(map_path,gene_name);
    /*output file path*/
    char *output_data_path=(char*)calloc(256,sizeof(char));
    output_data_path=strcat(output_data_path,output_path);

    beta_size_cnt(map_path,beta_size);
    nlex=beta_size[0];//THE COW OF BETA
    nclass=beta_size[1];//THE LINE OF BETA
    lda_nclass=nclass+1;//LDA模型计算添加一个主题。

    int emmax      = EMMAX_DEFAULT;		// default in lda.h
	int demmax     = DEMMAX_DEFAULT;	// default in lda.h
	double epsilon = EPSILON_DEFAULT;	// default in lda.h

	/* open data */
	if ((data = feature_matrix(gene_norm_data_path, &nlex, &dlenmax)) == NULL) {
		fprintf(stderr, "lda:: cannot open training data.\n");
		exit(1);
	}
	/* create the space for the sum of read*/
	 if ((sum_exon = (document*)calloc(1, sizeof(document))) == NULL)
    {
		fprintf(stderr,"gene_expression:: cannot allocate the total of read.\n");
		exit(1);
    }
	/* allocate parameters */
	if ((alpha = (double *)calloc(lda_nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate alpha.\n");
		exit(1);
	}
	/*LDA截取前面真实读段数据*/
	if ((short_alpha = (double *)calloc(nclass, sizeof(double))) == NULL) {
		fprintf(stderr, "lda:: cannot allocate alpha.\n");
		exit(1);
	}
	if ((beta = dmatrix(nlex, lda_nclass)) == NULL) {
		fprintf(stderr, "lda:: cannot allocate beta.\n");
		exit(1);
	}


	if((theta=(double *)calloc(nclass,sizeof(double)))==NULL)
	{
        fprintf(stderr, "lda:: cannot allocate theta.\n");
		exit(1);
	}
    if((isoform=dmatrix(nclass,nlex))==NULL)
    {
        fprintf(stderr, "lda:: cannot allocate isoform.\n");
		exit(1);
    }
    if((iso_expre=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
    {
        fprintf(stderr, "lda:: cannot allocate iso_expre.\n");
		exit(1);
    }
     if((iso_var=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
    {
        fprintf(stderr, "lda:: cannot allocate iso_var.\n");
		exit(1);
    }
    if((isolen=(double*)calloc(nclass,sizeof(double)))==NULL)//异构体长度
    {
        fprintf(stderr, "lda:: cannot allocate isolen.\n");
		exit(1);
    }

	read_beta(map_path,beta);
	normal_beta(beta,lda_nclass,nlex);


    /*lda learn */
	lda_learn (data, alpha, beta,lda_nclass, nlex, dlenmax,
		   emmax, demmax, epsilon);
    int cnt=0;
    double _sum=0;
    //alpha 归一化
    for(cnt=0;cnt<nclass;cnt++)
    {
        if(!isnan(alpha[cnt]))
        {
            short_alpha[cnt]=alpha[cnt];
        }
        else
        {
             short_alpha[cnt]=1.0/nclass;
        }

        _sum=short_alpha[cnt]+_sum;

    }
    for(cnt=0;cnt<nclass;cnt++)
    {
        short_alpha[cnt]=short_alpha[cnt]/_sum;
    }
    free_feature_matrix(data);
	/*calculate the gene expression*/
    data = feature_matrix(gene_data_path,&nlex, &dlenmax);
    sum_read_exon(data,sum_exon,nlex);//计算所有通道读段总数

    read_isolen(isoLen_path,isolen,nclass);

    double *iso_expre_sum;
    double *iso_expre_sum_2;
    if((iso_expre_sum=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
    {
        fprintf(stderr, "lda:: cannot allocate iso_expre.\n");
        exit(1);
    }
    if((iso_expre_sum_2=(double *)calloc(nclass,sizeof(double)))==NULL)//申请基因异构体表达值空间
    {
        fprintf(stderr, "lda:: cannot allocate iso_expre.\n");
        exit(1);
    }
    double gene_expre_sum=0;
    double gene_expre_sum_2=0;

    //采样，计算表达值
    const gsl_rng_type *T;
    gsl_rng *r ;
    gsl_rng_env_setup();
    T=gsl_rng_default;
    r=gsl_rng_alloc(T);
  /* end of GSL setup */


    for(cnt=0;cnt<SAMPLE_CNT;cnt++)
    {
        double t_gene_expre_sum=0;
        int t_cnt=0;
        gsl_ran_dirichlet(r,nclass,short_alpha,theta);
        //create_sampling(short_alpha,theta,nclass);//DIRCHLET sampling
        read_mapping(isoform,sum_exon,theta,beta,nclass,nlex);//读段映射
        expr_sig(isoform,nclass,nlex,isolen,iso_expre,SEQ_DEPTH);//
        //write_log("iso_expre_sample",iso_expre,nclass,"test_gene");
        for(t_cnt=0;t_cnt<nclass;t_cnt++)
        {
            t_gene_expre_sum=t_gene_expre_sum+iso_expre[t_cnt];//异构体表达值相加求基因表达值
            iso_expre_sum[t_cnt]=iso_expre_sum[t_cnt]+iso_expre[t_cnt];//每个异构体各自求和，计算方差用
            iso_expre_sum_2[t_cnt]=iso_expre_sum_2[t_cnt]+
                                        iso_expre[t_cnt]*iso_expre[t_cnt];//每个异构体和值
        }

        gene_expre_sum=gene_expre_sum+t_gene_expre_sum;
        gene_expre_sum_2=gene_expre_sum_2+gene_expre_sum*gene_expre_sum;
    }
    gsl_rng_free(r);
    for(cnt=0;cnt<nclass;cnt++)
    {
        iso_expre[cnt]=iso_expre_sum[cnt]/SAMPLE_CNT;
        iso_var[cnt]=fabs((iso_expre_sum_2[cnt]/SAMPLE_CNT)-(iso_expre[cnt]*iso_expre[cnt]));
    }

    double gene_expression;
    double gene_variance;

    gene_expression=gene_expre_sum/SAMPLE_CNT;
    gene_variance=fabs(gene_expre_sum_2/SAMPLE_CNT-(gene_expression*gene_expression));

    free(iso_expre_sum);
    free(iso_expre_sum_2);

   /* read_mapping(isoform,sum_exon,short_alpha,beta,nclass,nlex);//读段映射
    expr_sig(isoform,nclass,nlex,isolen,iso_expre,SEQ_DEPTH);
    for(cnt=0;cnt<nclass;cnt++)
    {
        iso_expre[cnt]=iso_expre[cnt]*1.0e9;
    }*/
    write_OUTCOME(output_data_path,iso_expre,iso_var,gene_expression,gene_variance,nclass,gene_name);

    free(map_path);
    free(gene_data_path);
    free(gene_norm_data_path);
    free(isoLen_path);
    free(output_data_path);

    free(iso_expre);
    free(iso_var);
    free(isolen);

	free_feature_matrix(data);
	free_sum_exon(sum_exon);
	free_dmatrix(beta, nlex);
    free(theta);
    free_dmatrix(isoform,nclass);
    free(alpha);
    free(short_alpha);
	//fclose(ap);
	//fclose(bp);
	return(0);
	//exit(0);
}



