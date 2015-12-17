/*
	J.Mater.Chem.,2008,18,4154¨C4164
	http://dx.doi.org/10.1039/B805147J 
*/

#pragma once
#include "Base.h"

struct TBMAtom{
	double A;
	double e;
	double p;
	double q;
	double r0;
};
typedef struct TBMAtom TBMATOM;

static TBMATOM TBMPt = {0.2975,2.695,10.612,4.004,2.7747},
			   TBMPd = {0.1746,1.718,10.867,3.742,2.7485},
			   TBMPtPd1 = {0.23,2.2,10.74,3.87,2.76},
			   TBMPtPd2 = {0.35,2.2,10.74,3.87,2.76},
			   TBMPtPd3 = {0.23,3.0,10.74,3.87,2.76};


double TBMCutEnergy(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N,double a0);

double TBMEnergy(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N);

TBMATOM GetTBMAtom(ATOM atom);
TBMATOM GetTwoTBMAtom(ATOM atom1,ATOM atom2);