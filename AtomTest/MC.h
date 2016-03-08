#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"

/***************************************
* @Name: ���ؿ��޷��㷨����
****************************************/
typedef struct
{
	int RSTEP;			//����������
	int temperature;	//�¶ȣ�������δʹ��
} MCPARA;

static const MCPARA _defMCPara = {
	100000,0
};

/***************************************
* @Name: ���ؿ��޷��㷨����
****************************************/
typedef struct
{
	int *note;			//ԭ������
	double E;			//����
} MCINDIVIDUAL;

/***************************************
* @Name: ���ؿ��޷��㷨ʵ��
****************************************/
typedef struct
{
	char *shape;		//��������	
	int N;				//ԭ������
	PE energyType;		//��������
	ALLOY alloy;		//�Ͻ�����
	ATOMNUM atomNum;	//���������ԭ�Ӹ���
	COOD cood;			//����ֵ
	COODDIS dis;		//ԭ�Ӽ����
	MCPARA para;		//�㷨����
	MCINDIVIDUAL one;	//����ʵ��	
	int clocks;			//��������
	int IJKL;			//�ҵ�����ֵ�Ĵ���
} MCInstance;


/***************************************
* @Name:MCPara_Init
* @Purpose����para��ֵĬ�ϵ�ֵ
			Ĭ��ֵΪ��̬���� _defMCPara
* @param��MCPARA *para --- ���ؿ����㷨����
* @return��void --- 
****************************************/
void MCPara_Init(MCPARA *para);

/***************************************
* @Name:MC2_InitWithMixing
* @Purpose�����ؿ��޶��Ͻ��ʼ�ṹΪ����ֲ�
* @param��char *shape --- ��״
            int N --- ��ԭ����
            int A --- ԭ��A����
            ATOM atomA --- ԭ��A����
            ATOM atomB --- ԭ��B����
			PE energy --- ����
            MCPARA *para --- ���ؿ��޷��㷨����
            char *Output --- �洢�ļ���
* @return��void --- 
****************************************/
void MC2_InitWithMixing(char *shape,int N,int A,
						ATOM atomA,ATOM atomB,
						PE energy,MCPARA *para,char *output);

/***************************************
* @Name:MC3_InitWithMixing
* @Purpose�����ؿ����㷨 ���Ͻ� ��ʼ�ṹΪ����ֲ�
* @param��char *shape --- ��״
            int N --- ��ԭ����
            int A --- ԭ��A����
            int B --- ԭ��B����
            ATOM atomA --- ԭ��A����
            ATOM atomB --- ԭ��B����
            ATOM atomC --- ԭ��C����
			PE energy --- ����
            MCPARA *para --- ���ؿ��޷��㷨����
            char *Output --- �洢�ļ�������
* @return��void --- 
****************************************/
void MC3_InitWithMixing(char *shape,int N,
						int A,int B,
						ATOM atomA,ATOM atomB,ATOM atomC,
						PE energy,MCPARA *para,char *output);

/***************************************
* @Name:MC_InitWithMixing
* @Purpose�����ؿ����㷨 ��Ͻ� ��ʼ�ṹ����Ų�
* @param��char *shape --- ��������
            int N --- ԭ���ܸ���
            ATOMNUM *atomNum --- �������ԭ������ 
            ALLOY *alloy --- �Ͻ�����
            PE energy --- ��������
            MCPARA *para --- ���ؿ����㷨����
            char *output --- �洢�ļ�������
* @return��void --- 
****************************************/
void MC_InitWithMixing(char *shape, int N, 
					   ATOMNUM *atomNum, ALLOY *alloy, 
					   PE energy, MCPARA *para, char *output);

/***************************************
* @Name:MC_Start
* @Purpose���������ؿ����㷨ʵ��������[output]�ļ��������ɽ���ļ��ͷ����ļ�
* @param��MCInstance *instance --- ���ؿ����㷨ʵ��
            char *output --- �洢�ļ�������
* @return��void --- 
****************************************/
void MC_Start(MCInstance *instance, char *output);


/***************************************
* @Name:MC_PrintMsg
* @Purpose����ӡ���ؿ����㷨ʵ������Ϣ
* @param��MCInstance *instance --- ���ؿ����㷨ʵ��
* @return��void --- 
****************************************/
void MC_PrintMsg( MCInstance *instance, char *output );

/***************************************
* @Name:MC_EnergyFile
* @Purpose����������ļ�����¼ÿһ����Ч������ֵ
			�ļ�·��Ϊ��[output]\\energy.txt
			���������֮����Ҫ�ֶ�fclose(fp)
* @param��MCInstance *instance --- ���ؿ����㷨ʵ��
            FILE **fp --- �ļ�ָ��ĵ�ַ
            char *output --- �ļ���·��
* @return��void --- 
****************************************/
void MC_EnergyFile(MCInstance *instance, FILE **fp, char *output);

/***************************************
* @Name:MC_ResultFile
* @Purpose�������ļ���·�������ɽ���ļ�:
			����ļ� - [output]\\result.txt
			��ͼ�ļ� - [output]\\Diamond.txt
* @param��MCInstance *instance --- ���ؿ����㷨ʵ��
            char *output --- �ļ���·��
* @return��void --- 
****************************************/
void MC_ResultFile(MCInstance *instance, char *output);

/***************************************
* @Name:MCIndividual_Init
* @Purpose����ʼ������ʵ����Ϊ�俪������Ŀռ�
			ÿһ������ʵ������Ҫ��ʼ��
* @param��MCINDIVIDUAL *one --- ���ؿ����㷨����ʵ��
            int N --- ��ԭ����
* @return��void --- 
****************************************/
void MCIndividual_Init(MCINDIVIDUAL *one,int N);

/***************************************
* @Name:MCIndividual_Free
* @Purpose���ͷŸ���ʵ����ʼ���п��ٵĿռ�
			��һ������ʵ��one���� GAIndividual_Init(one) ��
			����ʹ��ʱ�����ִ�� GAIndividual_Free(one)
* @param��MCINDIVIDUAL *one --- ���ؿ����㷨����ʵ��
* @return��void --- 
****************************************/
void MCIndividual_Free(MCINDIVIDUAL *one);



//����δ����,����ʹ��
//��ʼ��λ����
void MC2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//��ʼL0�ṹ�����ܲ�����
void MC2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);

//��ʼ�˿ǣ��ǲ�Ĳ���ΪsurfaceLayer 0:�����
void MC2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
void MC3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
void MC3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
