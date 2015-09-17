
#include <stdio.h>
#include <stdlib.h>
#include"read_beta.h"
#define BUFSIZE 65535
int beta_size_cnt(char* file_name,int *beta_size)
{
    FILE* fp;
    int col_cnt=0;
    int line_cnt=0;
    char line[BUFSIZE];
    //int line_cnt=0;
    if((fp=fopen(file_name,"r"))==NULL)
        return 0;
    while(!feof(fp))//read data until the end of the document
    {
        if((fgets(line,sizeof(line),fp))!=NULL)
        {
            int row=0;
            int cnt=0;
            while((line[row])!='\0')//read a line of data
            {
                if((line[row]!=0x0a)&&(line[row]!='\t'))//setp over the space
                {
                    cnt++;
                }
                row++;
            }
            line_cnt=cnt;
            col_cnt=col_cnt+1;
        }
    }
    fclose(fp);
    beta_size[0]=line_cnt;
    beta_size[1]=col_cnt;
    return(1);
}
int read_beta(char*file_name,double** Map)
{
    FILE* input;
    int theme_cnt=0;
    int word_cnt=0;
    char line[BUFSIZE];
    if((input=fopen(file_name,"r"))==NULL)
        return 0;
    while(!feof(input))//read data until the end of the document
    {
        if((fgets(line,sizeof(line),input))!=NULL)
        {
            int row=0;
            int cnt=0;
            while((line[row])!='\0')//read a line of data
            {
                if((line[row]!=0x0a)&&(line[row]!='\t'))//setp over the space
                {
                    Map[cnt][word_cnt]=line[row]-0x30;
                    //printf("%d\t",Map[cnt][word_cnt]);
                    cnt++;
                }
                row++;
            }
            theme_cnt=cnt;
        //printf("%s\n",line);
            word_cnt=word_cnt+1;
        }
    }
    fclose(input);
    int _cnt=0;
    for(_cnt=0;_cnt<theme_cnt;_cnt++)
    {
        Map[_cnt][word_cnt]=1;
    }
   //printf("theme count%d\n",theme_cnt);
    //printf("word count%d\n",word_cnt);
    return(1);
}
//归一化BETA矩阵，对每个主题归一化单词
int normal_beta(double** Map ,int nclass,int nlex)
{
    unsigned int cnt=0;
    for(cnt=0;cnt<nclass;cnt++)
    {
        double sum=0;
        unsigned int lex=0;
        for(lex=0;lex<nlex;lex++)
           // printf("%d\t",Map[cnt][lex]);
            sum=sum+Map[lex][cnt];
        for(lex=0;lex<nlex;lex++)
            Map[lex][cnt]=Map[lex][cnt]/sum;
    }
    return (1);
}
