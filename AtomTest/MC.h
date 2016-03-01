#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"


/***************************************
* @Name:MCS2_InitWithMixing
* @Purpose�����ؿ�����Ͻ��ʼ�ṹΪ����ֲ�
* @param��char *shape --- ��״
            int N --- ��ԭ����
            int A --- ԭ��A����
            ATOM atomA --- ԭ��A����
            ATOM atomB --- ԭ��B����
			PE energy --- ����
            int RSTEP --- ����������
            char *Output --- �洢�ļ���
* @return��void --- 
****************************************/
void MCS2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//��ʼ��λ����
void MCS2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//��ʼL0�ṹ�����ܲ�����
void MCS2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);


/***************************************
* @Name:MCS3_InitWithMixing
* @Purpose�����ؿ����㷨 ���Ͻ� ��ʼ�ṹΪ����ֲ�
* @param��char *shape --- ��״
            int N --- ��ԭ����
            int A --- ԭ��A����
            int B --- ԭ��B����
            ATOM atomA --- ԭ��A����
            ATOM atomB --- ԭ��B����
            ATOM atomC --- ԭ��C����
			PE energy --- ����
            int RSTEP --- ����������
            char *Output --- �洢�ļ�������
* @return��void --- 
****************************************/
void MCS3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
//��ʼ�˿ǣ��ǲ�Ĳ���ΪsurfaceLayer 0:�����
void MCS2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
void MCS3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
void MCS3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);

void MCS3_Start(char *shape,int N,int A,int B,int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
