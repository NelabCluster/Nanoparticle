#include "AEAM.h"

static double a[2] = {4.065,3.916};
static double n[2] = {1.41,1.43};
static double F0[2] = {2.969,4.403};
static double alpha[2] = {-4.968E-5, -7.447E-5};
static double k_1[2] = {54.351,98.633};
static double k0[2] = {-25.416,-46.656};
static double k1[2] = {-19.38,-31.402};
static double k2[2] = {-2.223,-13.393};
static double k3[2] = {0.573,10.527};
static double k4[2] = {-8.037,-17.908};
static double fe[2] = {3.87,5.64};
static double roue[2] = {49.3974,71.99};
static double Pe[2] = {180.3748,383.1};
static double kc = 0.75;
static double theta = 4.7;
static double mu = 1.14;
static double r1[2] = {2.8744,2.7690};
static double r2[2] = {4.0650,3.9160};
static double r3[2] = {4.9788,4.7963};
static double r4[2] = {5.7487,5.5380};
static double r5[2] = {6.4272,6.1916};
static double rce[2] = {6.2576,6.0282};
static double rc[2] = {5.5562,5.3526};
static double rab[2] = {2.9,2.834};
static double rc1 =  2.8217;
//	double r1[2] = {0.707*a[0],0.707*a[1]};
//	double r2[2] = {a[0],a[1]};
//	double r3[2] = {1.2248*a[0],1.2248*a[1]};
//	double r4[2] = {1.4142*a[0],1.4142*a[1]};
//	double r5[2] = {1.5811*a[0],1.5811*a[1]};
//	double rce[2] = {r4[0] + kc * (r5[0] - r4[0]),r4[1] + kc * (r5[1] - r4[1])};

double AEAMEnergy(int *note, double *R, ATOM atomA,ATOM atomB,ATOM atomC,int N)
{	
	int i,j;
	int atom,atom0,atom1;
	double r;
	double *f,*rou,*P,*F,*M;
	double Phi,E;
	f = calloc(N * N,sizeof(double));
	rou = calloc(N,sizeof(double));
	P = calloc(N,sizeof(double));
	F = calloc(N,sizeof(double));
	M = calloc(N,sizeof(double));
	for(i=0;i<N;i++)
	{
		for(j=0;j<N;j++)
		{	
			r = *(R+i*N+j);
			atom = *(note+i);
			if(i==j || r >rce[atom])
				*(f+i*N+j) = 0;
			else
				*(f+i*N+j) = fe[atom] * pow(r1[atom]/r,theta) * (rce[atom]-r)/(rce[atom]-r1[atom]) * (rce[atom]-r)/(rce[atom]-r1[atom]);
		}
	}
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
		{
			rou[i] += *(f+i*N+j);
			P[i] += *(f+i*N+j) * *(f+i*N+j);
		}

	for(i=0;i<N;i++) 
	{
		atom = note[i];
		F[i] = -F0[atom] * (1 - n[atom] * log(rou[i]/roue[atom])) * pow((rou[i] / roue[atom]),n[atom]);
		M[i] = alpha[atom] * (1-exp(-10000 * (log(P[i]/Pe[atom])) * (log(P[i]/Pe[atom]))));
	}

	Phi = 0;
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{	
			r = *(R+i*N+j);
			atom0 = note[i];
			atom1 = note[j];
			if(atom0 == atom1)
				Phi += PhiWithR(r,atom0);
			else
				Phi += 0.5 * mu * (PhiWithR(r*rab[atom0]/rc1,atom0) + PhiWithR(r*rab[atom1]/rc1,atom1)); 
		}
	}

	E = 0;
	for(i=0;i<N;i++)
	{
		E = E + F[i] + M[i];
	}
	E = E + Phi;

	free(f); 
	free(rou); 
	free(P); 
	free(F); 
	free(M);
	return E;
}

double PhiWithR(double r,int atom)
{
	double Phi;
	double temp;

	if(r > rc[atom])
		return 0;

	temp = 1 - r/r1[atom];
	Phi = k_1[atom] * exp(r1[atom]/r -1);
	Phi += k0[atom];
	Phi += k1[atom] * exp(temp);
	Phi += k2[atom] * exp(2 * temp);
	Phi += k3[atom] * exp(3 * temp);
	Phi += k4[atom] * exp(4 * temp);
	return Phi;
}

double AEAMCutEnergy(int *note, double *R, ATOM atom1,ATOM atom2,ATOM atom3,int N,double a0)
{
	return AEAMEnergy(note,R,atom1,atom2,atom3,N);
}