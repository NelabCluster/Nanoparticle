#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"


/***************************************
* @Name:MCS2_InitWithMixing
* @Purpose：蒙特卡洛二合金初始结构为随机分布
* @param：char *shape --- 形状
            int N --- 总原子数
            int A --- 原子A个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
			PE energy --- 势能
            int RSTEP --- 最大迭代次数
            char *Output --- 存储文件名
* @return：void --- 
****************************************/
void MCS2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//初始相位分离
void MCS2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//初始L0结构，可能不完整
void MCS2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);


/***************************************
* @Name:MCS3_InitWithMixing
* @Purpose：蒙特卡洛算法 三合金 初始结构为随机分布
* @param：char *shape --- 形状
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能
            int RSTEP --- 最大迭代次数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void MCS3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
//初始核壳，壳层的层数为surfaceLayer 0:最外层
void MCS2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
void MCS3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
void MCS3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);

void MCS3_Start(char *shape,int N,int A,int B,int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
