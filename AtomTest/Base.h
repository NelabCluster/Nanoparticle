#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include <direct.h>
#include <stdarg.h>

/***************************************
* @Purpose：用于可变函数中的结尾标识符
****************************************/
#define END -1

/***************************************
* @Purpose：产生[0,1)之间的实数
****************************************/
#define RAND1 rand()/(RAND_MAX+1.0)

/***************************************
* @Purpose：产生[0,N)之间的整数
****************************************/
#define RANDINT(N) (int)(N*rand()/(RAND_MAX+1.0))

/***************************************
* @Purpose：用于在DUBUG模式下打印日志
****************************************/
#ifdef _DUBUG
#define LOG(string) printf(s);
#else
#define LOG(string) 
#endif

/***************************************
* @Name:原子
****************************************/
enum Atom
{
	Ni,Cu,Rh,Pd,Ag,Ir,Pt,Au,Fe,Pd1,Pd2,Pd3
};
typedef enum Atom ATOM;

struct AtomPara
{
	double a;		//晶格常数
	char name[5];	//原子名称
};
typedef struct AtomPara ATOMPARA;

static ATOMPARA NiPara = {3.5157,"Ni"},
				CuPara = {3.6030,"Cu"},
				RhPara = {3.7984,"Rh"},
				PdPara = {3.8813,"Pd"},
				AgPara = {4.0691,"Ag"},
				IrPara = {3.8344,"Ir"},
				PtPara = {3.9163,"Pt"},
				AuPara = {4.0651,"Au"},
				FePara = {3.51,"Fe"},
				Pd1Para = {3.8813,"Pd"},
				Pd2Para = {3.8813,"Pd"},
				Pd3Para = {3.8813,"Pd"};

/***************************************
* @Name: 合金
****************************************/
struct Alloy
{
	ATOM *atoms;		//原子种类
	int atomTypeCount;	//原子种类数
};
typedef struct Alloy ALLOY;

typedef struct
{
	int *numberOfAtom;	//各种类原子数量
	int atomTypeCount;	//原子种类数
} ATOMNUM;

typedef struct 
{
	double *R;		//距离矩阵
	int N;			//原子数
	double Rmax;	//原子间的最大距离
	double Rmin;	//原子间的最小距离
} COODDIS;

typedef struct
{
	double *x;		//x轴坐标值
	double *y;		//y轴坐标值
	double *z;		//z轴坐标值
	int N;			//原子数量
} COOD;

/***************************************
* @Name:纳米粒子构型名称
****************************************/
static char SP[] = "SP",
			OH[] = "OH",
			TO[] = "TO",
			CU[] = "CU",
			RD[] = "RD",
			TH[] = "TH",
			DH[] = "DH",
			IH[] = "IH",
			THH210[] = "THH210",
			THH310[] = "THH310",
			THH520[] = "THH520",
			THH730[] = "THH730";


/***************************************
* @Name:势能
****************************************/
enum PotentialEnergy
{
	QSC,QSCNewModel,AEAM,Johnson,TBM
};

typedef enum PotentialEnergy PE;

typedef double (*PEnergy3) (int *note, double *R, ATOM atom1,ATOM atom2,ATOM atom3,int N);
typedef double (*PCutEnergy3) (int *note, double *R, ATOM atom1,ATOM atom2,ATOM atom3,int N,double a0);
typedef double (*PEnergy) (int *note, double *R, ATOM* atoms, int atomsN, int N);
typedef double (*PCutEnergy) (int *note, double *R, ATOM* atoms, int atomsN, int N);