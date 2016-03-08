#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include <direct.h>
#include <stdarg.h>

/***************************************
* @Purpose�����ڿɱ亯���еĽ�β��ʶ��
****************************************/
#define END -1

/***************************************
* @Purpose������[0,1)֮���ʵ��
****************************************/
#define RAND1 rand()/(RAND_MAX+1.0)

/***************************************
* @Purpose������[0,N)֮�������
****************************************/
#define RANDINT(N) (int)(N*rand()/(RAND_MAX+1.0))

/***************************************
* @Purpose��������DUBUGģʽ�´�ӡ��־
****************************************/
#ifdef _DUBUG
#define LOG(string) printf(s);
#else
#define LOG(string) 
#endif

/***************************************
* @Name:ԭ��
****************************************/
enum Atom
{
	Ni,Cu,Rh,Pd,Ag,Ir,Pt,Au,Fe,Pd1,Pd2,Pd3
};
typedef enum Atom ATOM;

struct AtomPara
{
	double a;		//������
	char name[5];	//ԭ������
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
* @Name: �Ͻ�
****************************************/
struct Alloy
{
	ATOM *atoms;		//ԭ������
	int atomTypeCount;	//ԭ��������
};
typedef struct Alloy ALLOY;

typedef struct
{
	int *numberOfAtom;	//������ԭ������
	int atomTypeCount;	//ԭ��������
} ATOMNUM;

typedef struct 
{
	double *R;		//�������
	int N;			//ԭ����
	double Rmax;	//ԭ�Ӽ��������
	double Rmin;	//ԭ�Ӽ����С����
} COODDIS;

typedef struct
{
	double *x;		//x������ֵ
	double *y;		//y������ֵ
	double *z;		//z������ֵ
	int N;			//ԭ������
} COOD;

/***************************************
* @Name:�������ӹ�������
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
* @Name:����
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