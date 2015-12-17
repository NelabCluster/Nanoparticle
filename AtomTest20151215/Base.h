#pragma once
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h> 
#include <string.h>
#include <direct.h>
#define RAND1 rand()/(RAND_MAX+1.0)

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
	double a;
	char name[5];
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
* @Name:纳米粒子形状
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