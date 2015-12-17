#include "QSC.h"

#define NN 100000
static double NUM1[NN],MUM1[NN];
static double NUM2[NN],MUM2[NN];
static double NUM3[NN],MUM3[NN];
static double NUM12[NN],MUM12[NN];
static double NUM13[NN],MUM13[NN];
static double NUM23[NN],MUM23[NN];
static double NUMP1[NN],MUMP1[NN];
static double NUMP2[NN],MUMP2[NN];
static double NUMP3[NN],MUMP3[NN];


void SetEnergyPow(ATOM atom1,ATOM atom2,ATOM atom3)
{
	int i;
	double XI;
	double n1,m1,n2,m2,n3,m3,n12,n13,n23,m12,m13,m23;
	QSCATOM atomA,atomB,atomC;

	atomA = GetQSCAtom(atom1); atomB = GetQSCAtom(atom2); atomC = GetQSCAtom(atom3);
	n1 = atomA.n; m1 = atomA.m; n2 = atomB.n; m2 = atomB.m;	n3 = atomC.n; m3 = atomC.m;
	n12 = (n1+n2)/2; m12 = (m1+m2)/2; n13 = (n1+n3)/2; m13 = (m1+m3)/2; n23 = (n2+n3)/2; m23 = (m2+m3)/2;

	// 将指数计算做成数组 
	for(i=0;i<NN;i++)
	{
		XI = (double)i/20000;
		NUM1[i] = pow(XI,n1); NUM2[i] = pow(XI,n2); NUM3[i] = pow(XI,n3);
		NUM12[i] = pow(XI,n12); NUM13[i] = pow(XI,n13); NUM23[i] = pow(XI,n23);
		MUM1[i] = pow(XI,m1); MUM2[i] = pow(XI,m2); MUM3[i] = pow(XI,m3);
		MUM12[i] = pow(XI,m12); MUM13[i] = pow(XI,m13); MUM23[i] = pow(XI,m23);
	}
}

//三合金能量
double QSCEnergy3(int *note,double *R,ATOM atom1,ATOM atom2,ATOM atom3,int N)
{
	int not,not1;
	int i,j;
	double R1;
	double RR,RR1,DRR1;
	int IRR1;
	double FMJV,FMJP;
	double EP=0;
    double *PEN,*VEN;
	double n1,m1,e1,c1,a1;
	double n2,m2,e2,c2,a2;
	double n3,m3,e3,c3,a3;
	double n12,n13,n23;
	double m12,m13,m23;
	double e12,e13,e23;
	double a12,a13,a23;
	QSCATOM atomA,atomB,atomC;

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	for(i=0;i<N;i++){
		PEN[i] = 0;
		VEN[i] = 0;
	}

	atomA = GetQSCAtom(atom1); atomB = GetQSCAtom(atom2); atomC = GetQSCAtom(atom3);
	n1 = atomA.n; m1 = atomA.m; e1 = atomA.e; a1 = atomA.a; c1 = atomA.c; 
	n2 = atomB.n; m2 = atomB.m; e2 = atomB.e; a2 = atomB.a; c2 = atomB.c; 
	n3 = atomC.n; m3 = atomC.m; e3 = atomC.e; a3 = atomC.a; c3 = atomC.c; 
	n12 = (n1+n2)/2; n13 = (n1+n3)/2; n23 = (n2+n3)/2;
	m12 = (m1+m2)/2; m13 = (m1+m3)/2; m23 = (m2+m3)/2;
	e12 = sqrt(e1*e2); e13 = sqrt(e1*e3); e23 = sqrt(e2*e3);
	a12 = (a1+a2)/2; a13 = (a1+a3)/2; a23 = (a2+a3)/2;
	
	for(i=0;i<N-1;i++){	
		not = *(note+i);
		for(j=i+1;j<N;j++)
		{   
			R1 = *(R+i*N+j);
            not1 = *(note+j);
			if(not1 == not)
			{
				if(not == 0)
				{
					RR = a1/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e1*(NUM1[IRR1]+(NUM1[IRR1+1]-NUM1[IRR1])*DRR1);
					FMJP = MUM1[IRR1]+(MUM1[IRR1+1]-MUM1[IRR1])*DRR1;
				} else if(not == 1)
				{
					RR = a2/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e2*(NUM2[IRR1]+(NUM2[IRR1+1]-NUM2[IRR1])*DRR1);
					FMJP = MUM2[IRR1]+(MUM2[IRR1+1]-MUM2[IRR1])*DRR1;
				} else{
					RR = a3/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e3*(NUM3[IRR1]+(NUM3[IRR1+1]-NUM3[IRR1])*DRR1);
					FMJP = MUM3[IRR1]+(MUM3[IRR1+1]-MUM3[IRR1])*DRR1;
				}
			}
			else
			{
				if((not+not1)==1)
				{
					
					RR = a12/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e12*(NUM12[IRR1]+(NUM12[IRR1+1]-NUM12[IRR1])*DRR1);
					FMJP = MUM12[IRR1]+(MUM12[IRR1+1]-MUM12[IRR1])*DRR1;
				} else if((not+not1)==2)
				{
					RR = a13/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e13*(NUM13[IRR1]+(NUM13[IRR1+1]-NUM13[IRR1])*DRR1);
					FMJP = MUM13[IRR1]+(MUM13[IRR1+1]-MUM13[IRR1])*DRR1;				
				} else
				{
					RR = a23/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e23*(NUM23[IRR1]+(NUM23[IRR1+1]-NUM23[IRR1])*DRR1);
					FMJP = MUM23[IRR1]+(MUM23[IRR1+1]-MUM23[IRR1])*DRR1;
				}
			}

		    VEN[i] = VEN[i]+FMJV*0.5;
			VEN[j] = VEN[j]+FMJV*0.5;
			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}
	}
	
	for(i=0;i<N;i++)
	{	
		not = *(note+i);
		if(not==0)
			PEN[i] = -c1*e1*sqrt(PEN[i]);
		else if(not==1)
			PEN[i] = -c2*e2*sqrt(PEN[i]);
		else
			PEN[i] = -c3*e3*sqrt(PEN[i]);
		
		EP = EP+VEN[i]+PEN[i];
	}
	

	free(PEN);
	free(VEN);
	return EP;
}

//截断的能量 三合金
double QSCCutEnergy3(int *note,double *R,ATOM atom1, ATOM atom2,ATOM atom3,int N,double a0)
{
	int not,not1;
	int i,j;
	double R1;
	double RR,RR1,DRR1;
	int IRR1;
	double FMJV,FMJP;
	double EP=0;
    double *PEN,*VEN;
	double n1,m1,e1,c1,a1;
	double n2,m2,e2,c2,a2;
	double n3,m3,e3,c3,a3;
	double n12,n13,n23;
	double m12,m13,m23;
	double e12,e13,e23;
	double a12,a13,a23;
	QSCATOM atomA,atomB,atomC;

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	
	atomA = GetQSCAtom(atom1); atomB = GetQSCAtom(atom2); atomC = GetQSCAtom(atom3);
	n1 = atomA.n; m1 = atomA.m; e1 = atomA.e; a1 = atomA.a; c1 = atomA.c; 
	n2 = atomB.n; m2 = atomB.m; e2 = atomB.e; a2 = atomB.a; c2 = atomB.c; 
	n3 = atomC.n; m3 = atomC.m; e3 = atomC.e; a3 = atomC.a; c3 = atomC.c; 
	n12 = (n1+n2)/2; n13 = (n1+n3)/2; n23 = (n2+n3)/2;
	m12 = (m1+m2)/2; m13 = (m1+m3)/2; m23 = (m2+m3)/2;
	e12 = sqrt(e1*e2); e13 = sqrt(e1*e3); e23 = sqrt(e2*e3);
	a12 = (a1+a2)/2; a13 = (a1+a3)/2; a23 = (a2+a3)/2;
	
	for(i=0;i<N-1;i++){	
		not = *(note+i);
		for(j=i+1;j<N;j++)
		{   
			R1 = *(R+i*N+j);
			if(R1 > 2*a0)
				continue;

            not1 = *(note+j);
			if(not1 == not)
			{
				if(not == 0)
				{
					RR = a1/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e1*(NUM1[IRR1]+(NUM1[IRR1+1]-NUM1[IRR1])*DRR1);
					FMJP = MUM1[IRR1]+(MUM1[IRR1+1]-MUM1[IRR1])*DRR1;
				} else if(not == 1)
				{
					RR = a2/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e2*(NUM2[IRR1]+(NUM2[IRR1+1]-NUM2[IRR1])*DRR1);
					FMJP = MUM2[IRR1]+(MUM2[IRR1+1]-MUM2[IRR1])*DRR1;
				} else{
					RR = a3/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e3*(NUM3[IRR1]+(NUM3[IRR1+1]-NUM3[IRR1])*DRR1);
					FMJP = MUM3[IRR1]+(MUM3[IRR1+1]-MUM3[IRR1])*DRR1;
				}
			}
			else
			{
				if((not+not1)==1)
				{
					RR = a12/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e12*(NUM12[IRR1]+(NUM12[IRR1+1]-NUM12[IRR1])*DRR1);
					FMJP = MUM12[IRR1]+(MUM12[IRR1+1]-MUM12[IRR1])*DRR1;
				} else if((not+not1)==2)
				{
					RR = a13/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e13*(NUM13[IRR1]+(NUM13[IRR1+1]-NUM13[IRR1])*DRR1);
					FMJP = MUM13[IRR1]+(MUM13[IRR1+1]-MUM13[IRR1])*DRR1;				
				} else
				{
					RR = a23/R1;
					RR1 = RR*20000;
					IRR1 = (int)(RR1);
					DRR1 = RR1-IRR1;
					FMJV = e23*(NUM23[IRR1]+(NUM23[IRR1+1]-NUM23[IRR1])*DRR1);
					FMJP = MUM23[IRR1]+(MUM23[IRR1+1]-MUM23[IRR1])*DRR1;
				}
			}

		    VEN[i] = VEN[i]+FMJV*0.5;
			VEN[j] = VEN[j]+FMJV*0.5;
			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}
	}

	for(i=0;i<N;i++)
	{	
		not = *(note+i);
		if(not==0)
			PEN[i] = -c1*e1*sqrt(PEN[i]);
		else if(not==1)
			PEN[i] = -c2*e2*sqrt(PEN[i]);
		else
			PEN[i] = -c3*e3*sqrt(PEN[i]);
		
		EP = EP+VEN[i]+PEN[i];
	}
	free(PEN);
	free(VEN);
	return EP;
}

QSCATOM GetQSCAtom(ATOM atom)
{
	QSCATOM zero = {0};
	switch(atom)
	{
	case Ni:
		return QSCNi;
	case Cu:
		return QSCCu;
	case Rh:
		return QSCRh;
	case Pd:
		return QSCPd;
	case Ag:
		return QSCAg;
	case Ir:
		return QSCIr;
	case Pt:
		return QSCPt;
	case Au:
		return QSCAu;
	}

	return zero;
}
