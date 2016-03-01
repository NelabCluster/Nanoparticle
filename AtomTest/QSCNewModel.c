#include "QSCNewModel.h"
#include "QSC.h"

static double *A00,*A10,*A11,*B00,*B10,*B01,*B11;

void setEnergyNewModel(ATOM atom1,ATOM atom2,double *R,int N)
{
	int i,j;
	double rij;
	double n1,m1,e1,c1,a1;
	double n2,m2,e2,c2,a2;
	double n3,m3,e3,a3;
	QSCATOM atomA,atomB;

	atomA = GetQSCAtom(atom1); atomB = GetQSCAtom(atom2);
	n1 = atomA.n; m1 = atomA.m; e1 = atomA.e; a1 = atomA.a; c1 = atomA.c; 
	n2 = atomB.n; m2 = atomB.m; e2 = atomB.e; a2 = atomB.a; c2 = atomB.c; 
	n3 = (n1 + n2) / 2; m3 = (m1 + m2) / 2; a3 = (a1 + a2)/2; e3 = sqrt(e1 * e2);
	

	A00 = calloc(N * N,sizeof(double));
	A10 = calloc(N * N,sizeof(double));
	A11 = calloc(N * N,sizeof(double));
	B00 = calloc(N * N,sizeof(double));
	B10 = calloc(N * N,sizeof(double));
	B01 = calloc(N * N,sizeof(double));
	B11 = calloc(N * N,sizeof(double));

	for(i=0;i<N-1;i++)
		for(j=i+1;j<N;j++)
		{
			rij = *(R + i*N + j);
			*(A00 + i*N +j) = e1 * pow(a1 / rij,n1); 
			*(A10 + i*N +j) = e3 * pow(a3 / rij,n3); 
			*(A11 + i*N +j) = e2 * pow(a2 / rij,n2); 
			*(B00 + i*N +j) = c1 * c1 * e1 * e1 * pow(a1 / rij,m1); 
			*(B10 + i*N +j) = c2 * c2 * e2 * e2 * pow(a3 / rij,m3);
			*(B01 + i*N +j) = c1 * c1 * e1 * e1 * pow(a3 / rij,m3);
			*(B11 + i*N +j) = c2 * c2 * e2 * e2 * pow(a2 / rij,m2);
		//	*(B00 + i*N +j) = pow(a1 / rij,m1); 
		//	*(B10 + i*N +j) = pow(a3 / rij,m3);
		//	*(B01 + i*N +j) = pow(a3 / rij,m3);
		//	*(B11 + i*N +j) = pow(a2 / rij,m2);
		}
}

void initModel()
{
	int i,j;
	int N = 443;
	char cood[30] ="THH210_443.txt";
	char num[100];
	FILE *fp;
	double *x,*y,*z,*R;
	double rij;
	double n1,m1,e1,c1,a1;
	double n2,m2,e2,c2,a2;
	double n3,m3,e3,a3;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	R = calloc(N*N,sizeof(double));
	
	//读坐标文件
	fp = fopen(cood,"r");
	fgets(num,99,fp);
	for(i=0;i<N;i++){
		fscanf(fp,"%lf %lf %lf",&x[i],&y[i],&z[i]);
	}
	fclose(fp);

	//计算距离
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);;
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
		}
	}

	//原子参数
	n1 = 11; m1 = 7; e1 = 9.7894E-3; a1 = 3.9163; c1 = 71.336; 
	n2 = 12; m2 = 6; e2 = 3.2864E-3; a2 = 3.8813; c2 = 148.205; 
	n3 = (n1 + n2) / 2; m3 = (m1 + m2) / 2; a3 = (a1 + a2)/2; e3 = sqrt(e1 * e2);
	

	A00 = calloc(N * N,sizeof(double));
	A10 = calloc(N * N,sizeof(double));
	A11 = calloc(N * N,sizeof(double));
	B00 = calloc(N * N,sizeof(double));
	B10 = calloc(N * N,sizeof(double));
	B01 = calloc(N * N,sizeof(double));
	B11 = calloc(N * N,sizeof(double));
	//计算系数矩阵
	for(i=0;i<N-1;i++)
		for(j=i+1;j<N;j++)
		{
			rij = *(R + i*N + j);
			*(A00 + i*N + j) = e1 * pow(a1 / rij,n1); 
			*(A10 + i*N + j) = e3 * pow(a3 / rij,n3); 
			*(A11 + i*N + j) = e2 * pow(a2 / rij,n2); 
			*(B00 + i*N + j) = pow(a1 / rij,m1); 
			*(B10 + i*N + j) = pow(a3 / rij,m3);
			*(B01 + i*N + j) = pow(a3 / rij,m3);
			*(B11 + i*N + j) = pow(a2 / rij,m2);

			*(A00 + j*N + i) = *(A00 + i*N +j); 
			*(A10 + j*N + i) = *(A10 + i*N +j); 
			*(A11 + j*N + i) = *(A11 + i*N +j); 
			*(B00 + j*N + i) = *(B00 + i*N +j); 
			*(B10 + j*N + i) = *(B10 + i*N +j);
			*(B01 + j*N + i) = *(B01 + i*N +j);
			*(B11 + j*N + i) = *(B11 + i*N +j);
		}

	free(x);
	free(y);
	free(z);
	free(R);
}

void freeSetEnergyNewModel()
{
	free(A00);
	free(A10);
	free(A11);
	free(B00);
	free(B10);
	free(B01);
	free(B11);
}

double QSCEnergyNewModel(int *note,ATOM atom1,ATOM atom2,int N)
{
	int i,j;
	int notei,notej;
	double E = 0;
	double *VEN,*PEN;
	double tempV,tempP;
	QSCATOM atomA,atomB;

	VEN = calloc(N,sizeof(double));
	PEN = calloc(N,sizeof(double));
	
	atomA = GetQSCAtom(atom1); atomB = GetQSCAtom(atom2);

	for(i=0;i<N-1;i++)
	{
		notei = note[i];
		for(j=i+1;j<N;j++)
		{
			notej = note[j];
			if(notei == 0 && notej == 0)
			{
				tempV = *(A00 + i* N + j);
				tempP = *(B00 + i* N + j);
			} else if(notei == 1 && notej == 0)
			{
				tempV = *(A10 + i* N + j);
				tempP = *(B10 + i* N + j);
			} else if(notei == 0 && notej == 1) 
			{
				tempV = *(A10 + i* N + j);
				tempP = *(B01 + i* N + j);
			} else
			{
				tempV = *(A11 + i* N + j);
				tempP = *(B11 + i* N + j);
			}

			VEN[i] += tempV * 0.5;
			VEN[j] += tempV * 0.5;
			PEN[i] += tempP;
			PEN[j] += tempP;

		}
	}

	

	for(i=0;i<N;i++)
	{
		if(note[i]==0)
			E = E + VEN[i]  - atomA.c * atomA.e * sqrt(PEN[i]);
		else
			E = E + VEN[i]  - atomB.c * atomB.e * sqrt(PEN[i]);
	}

	free(VEN);
	free(PEN);
	return E;
}


