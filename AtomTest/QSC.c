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

//三合金能量
double QSCEnergy(int *note, double *R, ATOM *atoms, int atomTypeCount, int N)
{
	int not,not1;
	int i,j;
	double R1;
	double RR;
	double FMJV,FMJP;
	double EP=0;
    double *PEN,*VEN;
	QSCATOM *QSCAtoms;

	QSCAtoms = calloc(atomTypeCount * atomTypeCount, sizeof(QSCATOM));
	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	
	for( i = 0; i < atomTypeCount; i++ )
	{
		for( j = 0; j < atomTypeCount; j++)
		{
			if( i == j)
				QSCAtoms[i*atomTypeCount + j] = GetQSCAtom(atoms[i]);
			else
				QSCAtoms[i*atomTypeCount + j] = GetQSCTwoAtom(atoms[i],atoms[j]);
		}
	}
	
	for(i=0;i<N-1;i++){	
		not = *(note+i);
		for(j=i+1;j<N;j++)
		{   
			R1 = *(R+i*N+j);
            not1 = *(note+j);

			RR = QSCAtoms[not * atomTypeCount + not1].a/R1;
			FMJV = QSCAtoms[not * atomTypeCount + not1].e*pow(RR,QSCAtoms[not * atomTypeCount + not1].n);
			FMJP = pow(RR,QSCAtoms[not * atomTypeCount + not1].m);

		    VEN[i] = VEN[i]+FMJV*0.5;
			VEN[j] = VEN[j]+FMJV*0.5;
			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}
	}
	
	for(i=0;i<N;i++)
	{	
		not = *(note+i);
		PEN[i] = - QSCAtoms[not * atomTypeCount + not].c
				 * QSCAtoms[not * atomTypeCount + not].e
				 * sqrt(PEN[i]);

		EP = EP+VEN[i]+PEN[i];
	}
	
	free(QSCAtoms);
	free(PEN);
	free(VEN);
	return EP;
}

//计算作用力
double QSCForce(int *note,double *R,
				double *x,double *y,double *z,
				double *FX,double *FY,double *FZ,
				ATOM *atoms, int atomTypeCount,int N)
{   
	int not1,not2;
	int i,j;
	double R2,RR;
	double FMJP,FMJ1,FMJ2,FK;
	double *PEN;
	double *FF;
	double FP=0;
	double maxF;
	QSCATOM *QSCAtoms,tempQSCAtom;

	QSCAtoms = calloc(atomTypeCount * atomTypeCount, sizeof(QSCATOM));
	PEN = calloc(N,sizeof(double));
	FF = calloc(N,sizeof(double));

	for( i = 0; i < atomTypeCount; i++ )
	{
		for( j = 0; j < atomTypeCount; j++)
		{
			if( i == j)
				QSCAtoms[i*atomTypeCount + j] = GetQSCAtom(atoms[i]);
			else
				QSCAtoms[i*atomTypeCount + j] = GetQSCTwoAtom(atoms[i],atoms[j]);
		}
	}

	maxF = 0;
	memset(FX,0,sizeof(double)*N);
	memset(FY,0,sizeof(double)*N);
	memset(FZ,0,sizeof(double)*N);

	for(i=0;i<N-1;i++)
	{	
		not1 = note[i];	
		for(j=i+1;j<N;j++)
		{   
			R2 = *(R+i*N+j); 
			not2 = note[j];
			
			RR = QSCAtoms[not1*atomTypeCount + not2].a/R2;
			FMJP = pow(RR,QSCAtoms[not1*atomTypeCount + not2].m);

			PEN[i] = PEN[i]+FMJP;
			PEN[j] = PEN[j]+FMJP;
		}
	}
	for(i=0;i<N;i++)
	{	
		not1 = note[i];
		PEN[i] = QSCAtoms[not1*atomTypeCount + not1].c*QSCAtoms[not1*atomTypeCount + not1].e
				 /2/sqrt(PEN[i]);
	}

	for(i=0;i<N-1;i++)
	{	
		not1 = note[i];	
		for(j=i+1;j<N;j++)
		{   
			R2 = *(R+i*N+j); 
			not2 = note[j];
			
			tempQSCAtom = QSCAtoms[not1*atomTypeCount + not2];
            RR = tempQSCAtom.a/R2;
			FMJ1 = pow(RR,tempQSCAtom.n+1);
			FMJ1 = -(tempQSCAtom.n*tempQSCAtom.e/tempQSCAtom.a)*FMJ1;
			FMJ2 = pow(RR,tempQSCAtom.m+1);
			FMJ2 = -(tempQSCAtom.m/tempQSCAtom.a)*FMJ2;
			
			FK = -(FMJ1 - (PEN[i] + PEN[j])*FMJ2);
			FX[i]=FX[i]+FK*(x[i]-x[j])/R2;
			FX[j]=FX[j]-FK*(x[i]-x[j])/R2;
			FY[i]=FY[i]+FK*(y[i]-y[j])/R2;
			FY[j]=FY[j]-FK*(y[i]-y[j])/R2;
			FZ[i]=FZ[i]+FK*(z[i]-z[j])/R2;
			FZ[j]=FZ[j]-FK*(z[i]-z[j])/R2;
		}
	}
	for(i=0;i<N;i++)
	{
		FF[i] = sqrt(FX[i]*FX[i]+FY[i]*FY[i]+FZ[i]*FZ[i]);
		FP = FP+FF[i]*FF[i];
		if(FF[i]>maxF)
			maxF = FF[i];
	}
	FP = sqrt(FP)/N;

	free(QSCAtoms);
	free(PEN);
	free(FF);
	return FP;
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

QSCATOM GetQSCTwoAtom(ATOM atom1, ATOM atom2)
{
	QSCATOM QSCAtom1, QSCAtom2;
	QSCATOM temp;

	QSCAtom1 = GetQSCAtom(atom1); QSCAtom2 = GetQSCAtom(atom2);
	temp.n = (QSCAtom1.n + QSCAtom2.n) / 2;
	temp.m = (QSCAtom1.m + QSCAtom2.m) / 2;
	temp.e = sqrt(QSCAtom1.e * QSCAtom2.e);
	temp.c = QSCAtom1.c;
	temp.a = (QSCAtom1.a + QSCAtom2.a) / 2;

	return temp;
}









