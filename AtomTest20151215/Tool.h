#pragma once
#include "Base.h"
      
double Rmax;
double Rmin;

ATOMPARA GetAtomPara(ATOM atom);

PEnergy3 GetEnergyFunction(PE type);
PCutEnergy3 GetCutEnergyFunction(PE type);


/***************************************
* @Name:ReadCood
* @Purpose：根据形状原子数读取坐标x,y,z
* @param：char *shape --- 形状 
            int N --- 总原子数
            double *x --- x轴坐标 指针，提前开辟空间
            double *y --- y轴坐标 指针，提前开辟空间
            double *z --- z轴坐标 指针，提前开辟空间
* @return：void --- 
****************************************/
void ReadCood(char *shape,int N,double *x,double *y,double *z);

/***************************************
* @Name:ReadFile
* @Purpose：根据文件读取原子分布note和坐标x,y,z
* @param：char *Input --- 文件名称
            int *note --- 原子分布 指针，提前开辟空间
            double *x --- x轴坐标 指针，提前开辟空间
            double *y --- y轴坐标 指针，提前开辟空间
            double *z --- z轴坐标 指针，提前开辟空间
            int N --- 总原子数
* @return：void --- 
****************************************/
void ReadFile(char *Input,int *note,double *x,double *y,double *z,int N);
//初始排列
void MixNoteInt3(int *note,int N,int A,int B);
void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//计算距离(x[],y[],z[],距离数组R,原子个数)
void Distance(double *x,double *y,double *z,double *R,int N);
//返回存储相对地址
char* StoragePath(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,char *Output);
//原子所在层
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);

void printDiamond3(int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,int N,char *path);
void printResult3(int *note,double *x,double *y,double *z,int N,double bili,char *path);
void printData(char *Line_Date,char *Line_End,int N);
//根据原子里质心距离排序坐标
int* orderCoodFromCore(double *x,double *y,double *z,int N);
int* OrderCoodAlongX(double *x,double *y,double *z,int N);
void L0NoteInt2(int *note,int N,int A,double *x,double *y,double *z);

/***************************************
* @Name:CoodNum3
* @Purpose：配位数上各原子的统计
* @param：char *input --- 输入文件
            int N --- 总原子数
            char *Output --- 文件存储路径
* @return：void --- 
****************************************/
void CoodNum3(char *input,int N,char *Output);


/***************************************
* @Name:pairCorrelation3
* @Purpose：径向分布函数
* @param：char *input --- 输入文件名称 
            int N --- 总原子数
            char *Output --- 文件存储路径
* @return：void --- 
****************************************/
void pairCorrelation3(char *input,int N,char *Output);

/***************************************
* @Name:ShellCount3
* @Purpose：各层原子统计
* @param：char *input --- 输入文件名称
            int N --- 总原子数
            char *Output --- 输出文件文件名称
* @return：void --- 
****************************************/
void ShellCount3(char *input,int N,char *Output);

//计算JA,JB,输入是二合金,结果打印出来.
void JAJB(char *input,int N);


/***************************************
* @Name:saveMatrix
* @Purpose：保存二维数组于本地
* @param：double *x --- 
            int N --- 
            char* output --- 
* @return：void --- 
****************************************/
void saveMatrix(double *x, int N, char* output);

double getLatticeParameter(ATOM atom1,ATOM atom2,ATOM atom3);

void SD_File(char *input, ATOM *atoms, int atomTypeCount, int N);