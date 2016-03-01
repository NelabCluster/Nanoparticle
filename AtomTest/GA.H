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
	int pc;						//交叉概率
	int pm;						//变异概率
	int convergenceGenerations;	//收敛代数
	char needOrderCood;			//是否需要排序；1：需要，0：不需要
	char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
	double adjustmentRate;		//调整概率，当needAdjustment是1时有效
} GAPARA;

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
* @Name:GA3_InitWithMixing
* @Purpose：遗传算法 三合金 初始结构随机排布
* @param：char *shape --- 形状
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能
            int POPSIZE --- 种群个体数量
            double pc --- 交叉概率
            double pm --- 变异概率
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);
void GA3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);

//初始核壳，壳层的层数为surfaceLayer 0:最外层
void GA3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);
void GA3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);

void GA3_Start(char *shape,int N,int A,int B,INDIVIDUAL *pop,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pm,double pc,double rate,char *Output);
void select_operator(struct individual *pop,int N,int POPSIZE,double bestEnergy,double worstEnergy);
void crossover_operator(struct individual *pop,int N,int POPSIZE,double pc);
void mutation_operator(struct individual *pop,int N,int POPSIZE,double pm);
void adjustment_operator(double rate, struct individual *pop, double *R, ATOM atomA, ATOM atomB, PE energy, ATOM atomC, int N, int POPSIZE);
void changeIndividual(struct individual *one,struct individual *two,int N);

//未使用
void GA_InitWithMixing(char *shape, int N, ATOMNUM atomNum, ALLOY alloy, PE energy, GAPARA para,char *output);
void GA_Start(char *shape, int N, ATOMNUM atomNum, INDIVIDUAL *pop, COOD cood, ALLOY alloy, PE energy, GAPARA para, char *output);