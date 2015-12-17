#pragma once
#include "Base.h"
/***************************************
* @Name:AEAMEnergy
* @Purpose：
* @param：int *note --- 
            double *R --- 
            int N --- 
* @return：double --- 
****************************************/
double AEAMEnergy(int *note, double *R, ATOM atom1,ATOM atom2,ATOM atom3,int N);
double PhiWithR(double r,int atom);

double AEAMCutEnergy(int *note, double *R, ATOM atom1,ATOM atom2,ATOM atom3,int N,double a0);//还是调用AEAMEnergy函数