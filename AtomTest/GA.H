#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"


/***************************************
* @Name: 遗传算法参数
****************************************/
typedef struct
{
	int popSize;				//种群个数
	double pc;					//交叉概率
	double pm;					//变异概率
	int convergenceGenerations;	//收敛代数
	char needOrderCood;			//是否需要排序；1：需要，0：不需要
	char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
	double adjustmentRate;		//调整概率，当needAdjustment是1时有效
} GAPARA;

/***************************************
* @Name: 遗传算法默认参数值
****************************************/
static const GAPARA _defGAPara = {
	100,0.8,0.0015,1000,1,1,1.0
};

/***************************************
* @Name: 遗传算法个体
****************************************/
struct individual
{	
	int number;			//保留，未使用
	int *chrom;			//原子排布序列		
	double energy;		//能量
	double value;		//适应值
};
typedef struct individual INDIVIDUAL;

/***************************************
* @Name: 遗传算法实例
****************************************/
typedef struct
{
	char *shape;		//构型
	int N;				//尺寸
	PE	energyType;		//势能类型
	ALLOY alloy;		//合金类型
	ATOMNUM atomNum;	//各类金属的原子个数
	COOD cood;			//坐标值
	COODDIS dis;		//原子间距离
	GAPARA para;		//遗传算法参数
	INDIVIDUAL *pop;	//种群
	INDIVIDUAL best;	//最优个体
	INDIVIDUAL worst;	//最差个体
} GAInstance;


/***************************************
* @Name:GA2_InitWithMixing
* @Purpose：遗传算法 二合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            PE energy --- 势能类型
            GAPARA *para --- 遗传算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,GAPARA *para,char *output);

/***************************************
* @Name:GA3_InitWithMixing
* @Purpose：遗传算法 三合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能类型
			GAPARA *para --- 遗传算法参数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,GAPARA *para,char *output);

//初始核壳，壳层的层数为surfaceLayer 0:最外层
void GA3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);
void GA3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);

void GA3_Start(char *shape,int N,int A,int B,INDIVIDUAL *pop,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pm,double pc,double rate,char *Output);
void select_operator(struct individual *pop,int N,int POPSIZE,double bestEnergy,double worstEnergy);
void crossover_operator(struct individual *pop,int N,int POPSIZE,double pc);
void mutation_operator(struct individual *pop,int N,int POPSIZE,double pm);
void adjustment_operator(double rate, struct individual *pop, double *R, ATOM atomA, ATOM atomB, ATOM atomC, PE energy,int N, int POPSIZE);
void adjustment_operator1(GAInstance *instance);

void GAIndividual_Init(INDIVIDUAL *one,int N);
void changeIndividual(struct individual *to,struct individual *from,int N);
void GAIndividual_Free(INDIVIDUAL *one);



void GAPara_Init(GAPARA *para);
void GA_InitWithMixing(char *shape, int N, ATOMNUM *atomNum, ALLOY *alloy, PE energy, GAPARA *para,char *output);
void GA_Start(GAInstance *instance, char *output);
void GA_FindBestAndWorst(GAInstance *instance);

void GA_PrintMsg(GAInstance *instance);

void GA_ResultFile(GAInstance *instance, char *output);