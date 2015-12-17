#include "TBM.h"

#define NNN 2

double TBMCutEnergy(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N,double a)
{
	int not;
	int i,j;
	double Rij;
	double FMJN,FMJV,FMJP;
	double EP=0;
    double *PEN,*VEN;
	double A1,e1,p1,q1,r01;
	double A2,e2,p2,q2,r02;
	double A3,e3,p3,q3,r03;
	TBMATOM atomA,atomB,atomAB;
	
	atomA = GetTBMAtom(atom1); atomB = GetTBMAtom(atom2); atomAB = GetTwoTBMAtom(atom1,atom2);
	A1 = atomA.A; e1 = atomA.e; p1 = atomA.p; q1 = atomA.q; r01 = atomA.r0;
	A2 = atomB.A; e2 = atomB.e; p2 = atomB.p; q2 = atomB.q; r02 = atomB.r0;
	A3 = atomAB.A; e3 = atomAB.e; p3 = atomAB.p; q3 = atomAB.q; r03 = atomAB.r0;
	

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	for(i=0;i<N-1;i++)
	{	
		not = *(note+i);	
		for(j=i+1;j<N;j++)
		{   
			Rij = *(R+i*N+j);
			if(Rij > NNN*a)
				continue;		
			if(*(note+j) == not)
			{
				if(not == 0)
				{
					FMJN = Rij/r01-1; 
					FMJV = A1*exp(-p1*FMJN);
					FMJP = e1*e1*exp(-2*q1*FMJN);
				}
				else
				{
					FMJN = Rij/r02-1;
					FMJV = A2*exp(-p2*FMJN);
					FMJP = e2*e2*exp(-2*q2*FMJN);
				}	
			}
			else
			{
				FMJN = Rij/r03-1;	
				FMJV = A3*exp(-p3*FMJN);
				FMJP = e3*e3*exp(-2*q3*FMJN);
			}

		    VEN[i] = VEN[i]+FMJV;
			VEN[j] = VEN[j]+FMJV;
			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}/* End of for(j=i+1;j<N;j++) */
	}/* End of for(i=0;i<N-1;i++) */

	for(i=0;i<N;i++)
	{	
		PEN[i] = sqrt(PEN[i]);		
		EP = EP+VEN[i]-PEN[i];
	}

	free(PEN);
	free(VEN);

	return EP;
}

/* 完整能量的计算 */
double TBMEnergy(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N)
{	
	int not;
	int i,j;
	double Rij;
	double FMJN,FMJV,FMJP;
	double EP=0;
    double *PEN,*VEN;
	double A1,e1,p1,q1,r01;
	double A2,e2,p2,q2,r02;
	double A3,e3,p3,q3,r03;
	TBMATOM atomA,atomB,atomAB;
	
	atomA = GetTBMAtom(atom1); atomB = GetTBMAtom(atom2); atomAB = GetTwoTBMAtom(atom1,atom2);
	A1 = atomA.A; e1 = atomA.e; p1 = atomA.p; q1 = atomA.q; r01 = atomA.r0;
	A2 = atomB.A; e2 = atomB.e; p2 = atomB.p; q2 = atomB.q; r02 = atomB.r0;
	A3 = atomAB.A; e3 = atomAB.e; p3 = atomAB.p; q3 = atomAB.q; r03 = atomAB.r0;

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));

	for(i=0;i<N-1;i++)
	{	
		not = *(note+i);	
		for(j=i+1;j<N;j++)
		{   
			Rij = *(R+i*N+j);		
			if(*(note+j) == not)
			{
				if(not == 0)
				{
					FMJN = Rij/r01-1; 
					FMJV = A1*exp(-p1*FMJN);
					FMJP = e1*e1*exp(-2*q1*FMJN);
				}
				else
				{
					FMJN = Rij/r02-1;
					FMJV = A2*exp(-p2*FMJN);
					FMJP = e2*e2*exp(-2*q2*FMJN);
				}	
			}
			else
			{
				FMJN = Rij/r03-1;	
				FMJV = A3*exp(-p3*FMJN);
				FMJP = e3*e3*exp(-2*q3*FMJN);
			}

		    VEN[i] = VEN[i]+FMJV;
			VEN[j] = VEN[j]+FMJV;
			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}/* End of for(j=i+1;j<N;j++) */
	}/* End of for(i=0;i<N-1;i++) */

	for(i=0;i<N;i++)
	{	
		PEN[i] = sqrt(PEN[i]);		
		EP = EP+VEN[i]-PEN[i];
	}

	free(PEN);
	free(VEN);

	return EP;
}

TBMATOM GetTBMAtom(ATOM atom)
{
	TBMATOM zero = {0};
	switch(atom)
	{
	case Pt:
		return TBMPt;
	case Pd1:
		return TBMPd;
	case Pd2:
		return TBMPd;
	case Pd3:
		return TBMPd;
	default:
		break;
	}
	return zero;
}

TBMATOM GetTwoTBMAtom(ATOM atom1,ATOM atom2)
{
	TBMATOM zero = {0};
	if(atom1 != Pt)
		return zero;
	switch(atom2)
	{
	case Pd1:
		return TBMPtPd1;
	case Pd2:
		return TBMPtPd2;
	case Pd3:
		return TBMPtPd3;
	default:
		break;
	}
	return zero;
}