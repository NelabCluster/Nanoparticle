#include "Johnson.h"

#define NN 100000
static double JSFeFr[NN],JSCuFr[NN];
static double JSFeTr[NN],JSCuTr[NN];
//double JSFeCuTr[NN];


void setupJohnson()
{
	int i;
	double XI;
	double FMJV1,FMJV2;

	//FeCu
/*	double re[2] = {2.481987,2.556162};
	double fe[2] = {1.885957,1.554485};
	double rhoe[2] = {20.041463,22.150141};
	double alpha[2] = {9.81870,7.669911};
	double beta[2] = {5.236411,4.090619};
	double A[2] = {0.392811,0.327584}; 
	double B[2] = {0.64624,0.468735};
	double kappa[2] = {0.170306,0.431307};
	double lambda[2] = {0.340613,0.86214};
	double Fn0[2] = {-2.534992,-2.176490};
	double Fn1[2] = {-0.059605,-0.140035};
	double Fn2[2] = {0.193065,0.285621};
	double Fn3[2] = {-2.282322,-1.750834};
	double F0[2] = {-2.54,-2.19};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,0.702991};
	double F3[2] = {-0.148770,0.683705};
	double eta[2] = {0.391750,0.921150};
	double Fe[2] = {-2.539945,-2.191675};*/
	//FePt
	double re[2] = {2.481987,2.771916};
	double fe[2] = {1.885957,2.336509};
	double rhoe[2] = {20.041463,34.108882};
	double alpha[2] = {9.81870,7.079952};
	double beta[2] = {5.236411,3.775974};
	double A[2] = {0.392811,0.449644}; 
	double B[2] = {0.64624,0.593713};
	double kappa[2] = {0.170306,0.413484};
	double lambda[2] = {0.340613,0.826967};
	double Fn0[2] = {-2.534992,-4.099542};
	double Fn1[2] = {-0.059605,-0.754764};
	double Fn2[2] = {0.193065,1.766503};
	double Fn3[2] = {-2.282322,-1.578274};
	double F0[2] = {-2.54,-4.17};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,3.474733};
	double F3[2] = {-0.148770,2.288323};
	double eta[2] = {0.391750,1.393490};
	double Fe[2] = {-2.539945,-4.174332};

	for(i=0;i<NN;i++)
	{	
		XI = (double)i/20000;
		FMJV1 = exp(- alpha[0] * (XI - 1)) / (1 + pow((XI - kappa[0]),20));
		FMJV2 = exp(- beta[0] * (XI - 1)) / (1 + pow((XI - lambda[0]),20));
		JSFeTr[i] = A[0] * FMJV1 - B[0] * FMJV2;
		JSFeFr[i] = fe[0] * FMJV2;

		FMJV1 = exp(- alpha[1] * (XI - 1)) / (1 + pow((XI - kappa[1]),20));
		FMJV2 = exp(- beta[1] * (XI - 1)) / (1 + pow((XI - lambda[1]),20));
		JSCuTr[i] = A[1] * FMJV1 - B[1] * FMJV2;
		JSCuFr[i] = fe[1] * FMJV2;

	//	JSFeCuTr[i] = (JSFeFr[i]/JSCuFr[i] * JSCuTr[i] + JSCuFr[i]/JSFeFr[i] * JSFeTr[i]) * 0.5;
	}	
}

double JohnsonCutEnergy(int *note,double *R,ATOM atomA,ATOM atomB,ATOM atomC,int N,double a0)
{
	int i,j;
	double r,rr,fa,fb;
	int notei,notej,not;
	double FMJV;
	double FMJPi,FMJPj;
	double FMJVa,FMJVb;
	double *PEN,*VEN;
	double rhon,rhoo,rho;
	double temp;
	double EP;
	double rr1;
	int indexrr;
	double DRR1;

	//FeCu
/*	double re[2] = {2.481987,2.556162};
	double fe[2] = {1.885957,1.554485};
	double rhoe[2] = {20.041463,22.150141};
	double alpha[2] = {9.81870,7.669911};
	double beta[2] = {5.236411,4.090619};
	double A[2] = {0.392811,0.327584}; 
	double B[2] = {0.64624,0.468735};
	double kappa[2] = {0.170306,0.431307};
	double lambda[2] = {0.340613,0.86214};
	double Fn0[2] = {-2.534992,-2.176490};
	double Fn1[2] = {-0.059605,-0.140035};
	double Fn2[2] = {0.193065,0.285621};
	double Fn3[2] = {-2.282322,-1.750834};
	double F0[2] = {-2.54,-2.19};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,0.702991};
	double F3[2] = {-0.148770,0.683705};
	double eta[2] = {0.391750,0.921150};
	double Fe[2] = {-2.539945,-2.191675};*/

	//FePt
	double re[2] = {2.481987,2.771916};
	double fe[2] = {1.885957,2.336509};
	double rhoe[2] = {20.041463,34.108882};
	double alpha[2] = {9.81870,7.079952};
	double beta[2] = {5.236411,3.775974};
	double A[2] = {0.392811,0.449644}; 
	double B[2] = {0.64624,0.593713};
	double kappa[2] = {0.170306,0.413484};
	double lambda[2] = {0.340613,0.826967};
	double Fn0[2] = {-2.534992,-4.099542};
	double Fn1[2] = {-0.059605,-0.754764};
	double Fn2[2] = {0.193065,1.766503};
	double Fn3[2] = {-2.282322,-1.578274};
	double F0[2] = {-2.54,-4.17};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,3.474733};
	double F3[2] = {-0.148770,2.288323};
	double eta[2] = {0.391750,1.393490};
	double Fe[2] = {-2.539945,-4.174332};

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	EP = 0;
	
	for(i=0;i<N-1;i++){	
		notei = *(note+i);
		for(j=i+1;j<N;j++)
		{   
			r = *(R+i*N+j);
            notej = *(note+j);
			if(r >= 3 * re[0])
				continue;
			if(notei == notej)
			{
				not = notei;
				rr = r / re[not];
				rr1 = rr * 20000;
				indexrr = (int)rr1;
				DRR1 = rr1 - indexrr;
				if(notei == 0){
					FMJV = JSFeTr[indexrr] + (JSFeTr[indexrr + 1]-JSFeTr[indexrr]) * DRR1;
					FMJPi = JSFeFr[indexrr] + (JSFeFr[indexrr + 1]-JSFeFr[indexrr]) * DRR1;
				}
				else {
					FMJV = JSCuTr[indexrr] + (JSCuTr[indexrr + 1]-JSCuTr[indexrr]) * DRR1;
					FMJPi = JSCuFr[indexrr] + (JSCuFr[indexrr + 1]-JSCuFr[indexrr]) * DRR1;
				}

				FMJPj = FMJPi;
			}
			else
			{
				rr = r / re[0];
				rr1 = rr * 20000;
				indexrr = (int)rr1;
				DRR1 = rr1 - indexrr;

				fa = JSFeFr[indexrr] + (JSFeFr[indexrr + 1]-JSFeFr[indexrr]) * DRR1;
				FMJVa = JSFeTr[indexrr] + (JSFeTr[indexrr + 1]-JSFeTr[indexrr]) * DRR1;
				
				rr = r / re[1];
				rr1 = rr * 20000;
				indexrr = (int)rr1;
				DRR1 = rr1 - indexrr;

				fb = JSCuFr[indexrr] + (JSCuFr[indexrr + 1]-JSCuFr[indexrr]) * DRR1;
				FMJVb = JSCuTr[indexrr] + (JSCuTr[indexrr + 1]-JSCuTr[indexrr]) * DRR1;

				FMJV = fb / fa * FMJVa + fa / fb * FMJVb;
				FMJV = FMJV / 2;

				FMJPi = (notei == 0)?fa:fb;
				FMJPj = (notej == 0)?fa:fb;
			}

		    VEN[i] = VEN[i]+FMJV*0.5;
			VEN[j] = VEN[j]+FMJV*0.5;
			PEN[i] = PEN[i]+FMJPj;
			PEN[j] = PEN[j]+FMJPi;
		}
	}

	for(i=0;i<N;i++)
	{	
		not = *(note+i);
		rhon = 0.85 * rhoe[not];
		rhoo = 1.15 * rhoe[not];
		
		
		rho = PEN[i];
		if(rho < rhon)
		{
			temp  = rho / rhon - 1;
			PEN[i] = Fn0[not] + Fn1[not] * temp;
			PEN[i] += Fn2[not] * temp * temp;
			PEN[i] += Fn3[not] * temp * temp * temp;
		} else if(rho < rhoo)
		{
			temp  = rho / rhoe[not] - 1;
			PEN[i] = F0[not] + F1[not] * temp;
			PEN[i] += F2[not] * temp * temp;
			PEN[i] += F3[not] * temp * temp * temp;
		} else 
		{
			temp  = rho / rhoe[not] - 1;
			PEN[i] = Fe[not] * (1 - eta[not]*log(temp)) * pow(temp, eta[not]);
		}
		
		EP = EP+VEN[i]+PEN[i];
	}
	free(PEN);
	free(VEN);
	return EP;
}

//二合金能量
double JohnsonEnergy(int *note,double *R,ATOM atomA,ATOM atomB,ATOM atomC,int N)
{
	int i,j;
	double r,rr,ff,fa,fb;
	int notei,notej,not;
	double FMJV1,FMJV2,FMJV;
	double FMJPi,FMJPj;
	double FMJVa1,FMJVa2,FMJVa,FMJVb1,FMJVb2,FMJVb;
	double *PEN,*VEN;
	double rhon,rhoo,rho;
	double temp;
	double EP;

	//FeCu
/*	double re[2] = {2.481987,2.556162};
	double fe[2] = {1.885957,1.554485};
	double rhoe[2] = {20.041463,22.150141};
	double alpha[2] = {9.81870,7.669911};
	double beta[2] = {5.236411,4.090619};
	double A[2] = {0.392811,0.327584}; 
	double B[2] = {0.64624,0.468735};
	double kappa[2] = {0.170306,0.431307};
	double lambda[2] = {0.340613,0.86214};
	double Fn0[2] = {-2.534992,-2.176490};
	double Fn1[2] = {-0.059605,-0.140035};
	double Fn2[2] = {0.193065,0.285621};
	double Fn3[2] = {-2.282322,-1.750834};
	double F0[2] = {-2.54,-2.19};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,0.702991};
	double F3[2] = {-0.148770,0.683705};
	double eta[2] = {0.391750,0.921150};
	double Fe[2] = {-2.539945,-2.191675};*/

	//FePt
	double re[2] = {2.481987,2.771916};
	double fe[2] = {1.885957,2.336509};
	double rhoe[2] = {20.041463,34.108882};
	double alpha[2] = {9.81870,7.079952};
	double beta[2] = {5.236411,3.775974};
	double A[2] = {0.392811,0.449644}; 
	double B[2] = {0.64624,0.593713};
	double kappa[2] = {0.170306,0.413484};
	double lambda[2] = {0.340613,0.826967};
	double Fn0[2] = {-2.534992,-4.099542};
	double Fn1[2] = {-0.059605,-0.754764};
	double Fn2[2] = {0.193065,1.766503};
	double Fn3[2] = {-2.282322,-1.578274};
	double F0[2] = {-2.54,-4.17};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,3.474733};
	double F3[2] = {-0.148770,2.288323};
	double eta[2] = {0.391750,1.393490};
	double Fe[2] = {-2.539945,-4.174332};

	PEN = calloc(N,sizeof(double));
	VEN = calloc(N,sizeof(double));
	EP = 0;
	
	for(i=0;i<N-1;i++){	
		notei = *(note+i);
		for(j=i+1;j<N;j++)
		{   
			r = *(R+i*N+j);
            notej = *(note+j);

			if(notei == notej)
			{
				not = notei;
				rr = r / re[not];
				FMJV1 = exp(- alpha[not] * (rr - 1)) / (1 + pow((rr - kappa[not]),20));
				FMJV2 = exp(- beta[notej] * (rr - 1)) / (1 + pow((rr - lambda[notej]),20));
				FMJV = A[not] * FMJV1 - B[not] * FMJV2;
				
				FMJPi = fe[not] * FMJV2;
				FMJPj = FMJPi;
			}
			else
			{
				rr = r / re[0];
				ff = exp(- beta[0] * (rr - 1)) / (1 + pow((rr - lambda[0]),20));
				fa = fe[0] * ff;
				FMJVa1 = exp(- alpha[0] * (rr - 1)) / (1 + pow((rr - kappa[0]),20));
				FMJVa2 = ff;
				FMJVa = A[0] * FMJVa1 - B[0] * FMJVa2;
				
				rr = r / re[1];
				ff = exp(- beta[1] * (rr - 1)) / (1 + pow((rr - lambda[1]),20));
				fb = fe[1] * ff;
				FMJVb1 = exp(- alpha[1] * (rr - 1)) / (1 + pow((rr - kappa[1]),20));
				FMJVb2 = ff;
				FMJVb = A[1] * FMJVb1 - B[1] * FMJVb2;

				FMJV = fb / fa * FMJVa + fa / fb * FMJVb;
				FMJV = FMJV / 2;

				FMJPi = (notei == 0)?fa:fb;
				FMJPj = (notej == 0)?fa:fb;
			}

		    VEN[i] = VEN[i]+FMJV*0.5;
			VEN[j] = VEN[j]+FMJV*0.5;
			PEN[i] = PEN[i]+FMJPj;
			PEN[j] = PEN[j]+FMJPi;
		}
	}

	for(i=0;i<N;i++)
	{	
		not = *(note+i);
		rhon = 0.85 * rhoe[not];
		rhoo = 1.15 * rhoe[not];
		
		
		rho = PEN[i];
		if(rho < rhon)
		{
			temp  = rho / rhon - 1;
			PEN[i] = Fn0[not] + Fn1[not] * temp;
			PEN[i] += Fn2[not] * temp * temp;
			PEN[i] += Fn3[not] * temp * temp * temp;
		} else if(rho < rhoo)
		{
			temp  = rho / rhoe[not] - 1;
			PEN[i] = F0[not] + F1[not] * temp;
			PEN[i] += F2[not] * temp * temp;
			PEN[i] += F3[not] * temp * temp * temp;
		} else 
		{
			temp  = rho / rhoe[not] - 1;
			PEN[i] = Fe[not] * (1 - eta[not]*log(temp)) * pow(temp, eta[not]);
		}
		
		EP = EP+VEN[i]+PEN[i];
	}
	free(PEN);
	free(VEN);
	return EP;
}