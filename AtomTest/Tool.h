#pragma once
#include "Base.h"
      
double Rmax;
double Rmin;

ATOMPARA GetAtomPara(ATOM atom);

PEnergy3 GetEnergyFunction(PE type);
PCutEnergy3 GetCutEnergyFunction(PE type);
PEnergy GetEnergyFunction1(PE type);

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
void ReadCood1(char *shape,int N,COOD *cood);

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
void MixNoteInt(int *note, int N, ATOMNUM *atomNum);

void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//计算距离(x[],y[],z[],距离数组R,原子个数)
void Distance(double *x,double *y,double *z,double *R,int N);
void Distance1(COOD *cood,COODDIS *dis);
//返回存储相对地址
char* StoragePath(char *shape,int N,ALLOY *alloy,ATOMNUM *atomNum,char *Output);
//原子所在层
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);

void printResult3(int *note,double *x,double *y,double *z,int N,double bili,char *path);
void printResult(int *note,int N,COOD *cood,char *path);
void printDiamond3(int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,int N,char *path);
void printDiamond(int *note,int N,COOD *cood,ALLOY *alloy,char *path);

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
            char *Output --- 输出文件名称
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

/***************************************
* @Name:getLatticeParameter3
* @Purpose：计算三合金的晶格
			如果计算二合金，只要是atom2 = atom3 即可
* @param：ATOM atom1 --- 原子类型1
            ATOM atom2 --- 原子类型2
            ATOM atom3 --- 原子类型3
* @return：double --- 晶格常数
****************************************/
double getLatticeParameter3(ATOM atom1,ATOM atom2,ATOM atom3);

/***************************************
* @Name:getLatticeParameter
* @Purpose：计算原子类型大于三的合金晶格
* @param：ATOM* atoms --- 原子类型
            int atomTypeCount --- 原子
* @return：double --- 晶格
****************************************/
double getLatticeParameter(ALLOY *alloy);

void SD_File(char *input, ATOM *atoms, int atomTypeCount, int N);

void Alloy_Init(ALLOY *alloy,...);
void Alloy_Copy(ALLOY *to, ALLOY *from);
void Alloy_Free(ALLOY *alloy);
void AtomNum_Init(ATOMNUM *atomNum,...);
void AtomNum_Copy(ATOMNUM *to, ATOMNUM *from);
void AtomNum_Free(ATOMNUM *atomNum);

void Cood_Init(COOD *cood,int N);
void Cood_Free(COOD *cood);