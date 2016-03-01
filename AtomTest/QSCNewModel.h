#pragma once
#include "Base.h"

/***************************************
* @Name:setEnergyNewModel
* @Purpose£º
* @param£ºATOM atomA --- 
            ATOM atomB --- 
            double *R --- 
            int N --- 
* @return£ºvoid --- 
****************************************/
void setEnergyNewModel(ATOM atom1,ATOM atom2,double *R,int N);
double QSCEnergyNewModel(int *note,ATOM atom1,ATOM atom2,int N);
void freeSetEnergyNewModel();
void initModel();