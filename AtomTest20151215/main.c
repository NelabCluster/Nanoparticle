#pragma once
/**************************势能原子***********************************/
//	QSC: Ni Cu Rh Pd Ag Ir Pt Au  
//	QSCNewModel: 占时不用
//	AEAM:
//	Johnson: Fe Pt Cu
//	TBM(Gupta): Pt Pd1 Pd2 Pd3
/**************************形状尺寸***********************************/
//	SP:	19 141 459 1061 2171 3781 5895 8685 12479	
//	CU: 63 365 1099 2457 4631 7813 12195
//	RD: 19 93 279 617 1163 1957 3055 4497 6339 8621 11399
//	TH: 35 165 455 969 1771 2925 4495 6545 9139 12341
//	OH: 19 85 231 489 891 1469 2255 3281 4579 6181 8119 10425
//	TO: 63 309 711 1289 2075 3101 4399 6001 7939 10245
//	DH: 104 608 1832 4096 7720 13024
//	IH: 308 2056 6524 14992
//	THH210: 63 443 1417 3285 6323 10831 25385 49251
//	THH310: 1249 9553 31825 74977
#include "Base.h"		//基本数据的定义
#include "Tool.h"		//通用函数的定义
#include "Algorithm.h"	//算法的定义
#include "Energy.h"		//势能的定义

void test()
{
	double *x,*y,*z,*R;
	int *note;
	int N;
	double E1,E2,E3;
	ATOM atoms[3] = {Pt,Pd,Rh};

	N = 1417;
	note = calloc(N,sizeof(int));
	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	R = calloc(N*N,sizeof(double));
	ReadFile("result.txt",note,x,y,z,N);
	Distance(x,y,z,R,N);
	SetEnergyPow(Pt,Pd,Rh);
	E1 = QSCEnergy3(note,R,Pt,Pd,Rh,N);
	E2 = QSCEnergy(note,R,atoms,3,N);
	printf("%f\n%f\n%f\n",E1,E2,E3);
	system("pause");
}

void test1()
{
	ATOM atoms[3] = {Pt,Pd,Rh};
	int N = 1417;

	SD_File("result.txt",atoms,3,N);
}

main()
{ 
	test1();

	//算法：蒙特卡洛算法；
	//初始结构：随机
	//形状：立法体
	//总原子数: 365
	//第一种原子的个数(这里就是Au):180
	//合金类型: Au-Ag;
	//势能：QSC
	//代数：10000
	//存储的文件夹名称：“MCS_MIX”，可以自己设定
//	MCS2_InitWithMixing("CU",365,180,Au,Ag,QSC,100000,"MCS_MIX");

	/****其他初始结构*****/
	//初始结构：相位分离
	//存储的文件夹名称：“MCS_PHASE”
//	MCS2_InitWithPhase("CU",365,180,Au,Ag,QSC,10000,"MCS_PHASE");
	//初始结构：核壳结构，第一种原子在内部
	//存储的文件夹名称：“MCS_CORE”
//	MCS2_InitWithFromCore("CU",365,180,Au,Ag,QSC,10000,"MCS_CORE");
	//初始结构: 有序初始结构
	//存储的文件夹名称：“MCS_L0”
//	MCS2_InitWithL0("CU",365,180,Au,Ag,QSC,10000,"MCS_L0");


	//算法：遗传算法；
	//初始结构：随机
	//形状：立法体
	//总原子数: 365
	//第一种原子的个数(这里就是Au):180
	//合金类型: Au-Ag
	//势能：QSC
	//种群大小：100
	//交叉概率：0.8
	//变异率：0.0015
	//调整概率：1.0
	//存储的文件夹名称：“GA_MIX”
//	GA2_InitWithMixing("CU",365,180,Au,Ag,QSC,100,0.8,0.0015,1.0,"GA_MIX");

	//算法：粒子群算法；
	//初始结构：随机
	//形状：24面体
	//总原子数: 3285
	//第一种原子的个数(这里就是Pt):1642
	//合金类型: Pt-Pd
	//势能：TBM
	//种群大小：100
	//惯性系数：0.25
	//调整率：1.0
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd1,TBM,100,0.25,1.0,"PSO_PARA1_YES",PSOYesImprove);//改进后PSO，第一组参数
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd2,TBM,100,0.25,1.0,"PSO_PARA2_YES",PSOYesImprove);//改进后PSO，第二组参数
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd3,TBM,100,0.25,1.0,"PSO_PARA3_YES",PSOYesImprove);//改进后PSO，第三组参数

//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd1,TBM,100,0.25,1.0,"PSO_PARA1_NO",PSONoImprove);//改进前PSO，第一组参数
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd2,TBM,100,0.25,1.0,"PSO_PARA2_NO",PSONoImprove);//改进前PSO，第二组参数
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd3,TBM,100,0.25,1.0,"PSO_PARA3_NO",PSONoImprove);//改进前PSO，第三组参数

}