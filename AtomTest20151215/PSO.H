#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"


enum PSOImprove
{
	PSONoImprove,
	PSOYesImprove
};
typedef enum PSOImprove PSOIMPROVE;


/***************************************
* @Name: 交换子
****************************************/
struct PSOSwap
{
	int element[2];			
	struct PSOSwap *next;	
};
typedef struct PSOSwap PSOSWAP;


/***************************************
* @Name:粒子群个体
****************************************/
struct PSOIndividual
{
	int *chrom;						//解
	double energy;					//能量
	PSOSWAP *swapList;				//交换序
	struct PSOIndividual *pbest;	//个体最优
};
typedef struct PSOIndividual PSOINDIVIDUAL;

/***************************************
* @Name:PSO2_InitWithMixing
* @Purpose：粒子群算法二合金初始结构随机分布
* @param：char *shape --- 形状
            int N --- 原子总数
            int A --- 原子A个数
            ATOM atomA --- 原子A
            ATOM atomB --- 原子B
			PE energy --- 势能
            int Popsize --- 个体数
            double w --- 惯性常数
            double rate --- 变异率
            char *Output --- 文件存储路径
* @return：void --- 
****************************************/
void PSO2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output,PSOIMPROVE improve); 
void PSO2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output);
void PSO2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output);
void PSO2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output);

void PSO2_Start(char *shape,int N,int A,PSOINDIVIDUAL *pop,double *x,double *y,double *z,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output,PSOIMPROVE improve);
void updateSwapList(double w, PSOINDIVIDUAL *pop, PSOINDIVIDUAL *gbest, int POPSIZE, int N);
void upDateChorm(PSOINDIVIDUAL *pop, int POPSIZE);
void PSOMutation(double rate,double *R, PSOINDIVIDUAL *pop, ATOM atomA, ATOM atomB, PE energy, int POPSIZE, int N);
void updatePbest(double *R, PSOINDIVIDUAL *pop, ATOM atomA, ATOM atomB, PE energy, int POPSIZE, int N);
void updateGbest(PSOINDIVIDUAL **best_pop,PSOINDIVIDUAL *pop,int POPSIZE);

PSOSWAP* initSwapSequence(int size,int N);
PSOSWAP* creatSwap(int element1,int element2);
void changePSOIndividual(PSOINDIVIDUAL *one,PSOINDIVIDUAL *two,int N);
PSOSWAP* swapSubtraction(PSOINDIVIDUAL *best, PSOINDIVIDUAL *pop, int N);
PSOSWAP* swapAddition(double w, PSOINDIVIDUAL *pop, PSOSWAP *swapList1, PSOSWAP *swapList2);
PSOSWAP* swapMultiplication(double rate, PSOSWAP *swapList, PSOSWAP *newSwapList);
void printfPSOSwap(PSOSWAP *head);
void printfPSOIndividual(PSOINDIVIDUAL *individual,int N);