/*
	The Quantum Sutton-Chen Many-Body Potential for Properties of fcc Metals 
	https://www.researchgate.net/publication/2327883_The_Quantum_Sutton-Chen_Many-Body_Potential_for_Properties_of_fcc_Metals
*/
#pragma once
#include "Base.h"

struct QSCAtom{
	double n;
	double m;
	double e;
	double c;
	double a;
};
typedef struct QSCAtom QSCATOM;

static QSCATOM QSCNi = {10,5,7.3767E-3,84.745,3.5157},
			   QSCCu = {10,5,5.7921E-3,84.843,3.6030},
			   QSCRh = {13,5,2.4612E-3,305.499,3.7984},
			   QSCPd = {12,6,3.2864E-3,148.205,3.8813},
			   QSCAg = {11,6,3.9450E-3,96.524,4.0691},
			   QSCIr = {13,6,3.7674E-3,224.815,3.8344},
			   QSCPt = {11,7,9.7894E-3,71.336,3.9163},
			   QSCAu = {11,8,7.8052E-3,53.581,4.0651};




//计算QSC势能，使用查表法的时候的初始化函数
void SetEnergyPow(ATOM atom1,ATOM atom2,ATOM atom3);
//计算QSC势能能量

double QSCEnergy3(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N);
double QSCCutEnergy3(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N,double a0);

double QSCEnergy(int *note, double *R, ATOM *atom, int atomTypeCount, int N);
double QSCCutEnergy(int *note, double *R, ATOM *atoms, int atomTypeCount, int N, double a0);


double QSCForce(int *note,double *R,
				double *x,double *y,double *z,
				double *FX,double *FY,double *FZ,
				ATOM *atoms, int atomTypeCount,int N);

//	double QSCCutForce();


QSCATOM GetQSCAtom(ATOM atom);
QSCATOM GetQSCTwoAtom(ATOM atom1, ATOM atom2);
