#pragma once
#include "Base.h"



void setupJohnson();

double JohnsonCutEnergy(int *note,double *R,ATOM atomA,ATOM atomB,ATOM atomC,int N,double a0);

double JohnsonEnergy(int *note,double *R,ATOM atomA,ATOM atomB,ATOM atomC,int N);

