#include "MC.h"

//蒙特卡洛算法 二合金
void MCS2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	if(A > N)
	{
		printf("原子个数超过总原子数量");
		return;
	}
	MCS3_InitWithMixing(shape,N,A,N-A,atomA,atomB,atomB,energy,RSTEP,Output);
}

void MCS2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	if(A > N)
	{
		printf("原子个数超过总原子数量");
		return;
	}
	MCS3_InitWithFromCore(shape,N,A,N-A,atomA,atomB,atomB,energy,RSTEP,Output);
}

void MCS2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	int *note;
	double *x,*y,*z;
	char *Line_End;
	char Line_Date[100];
	FILE *fp;

	if(A > N)
	{
		printf("原子个数超过总原子数量");
		return;
	}

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	Line_End = StoragePath(shape,N,A,N-A,atomA,atomB,atomB,Output);
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	//判断是否存在result.txt，存在则继续实验，不存在则到Cood文件夹读取坐标然后产生随机排布。
	if((fp=fopen(Line_Date,"r"))==NULL){
		printf("\nfile   not   exist\n");
		ReadCood(shape,N,x,y,z);
		srand((unsigned)time(NULL));
		PhaseSeparationNoteInt2(note,N,A,x,y,z);
	} else     
	{
		printf("\nfile   exist\n");  
		ReadFile(Line_Date,note,x,y,z,N);
		fclose(fp);
	}
	
	MCS3_Start(shape,N,A,N-A,note,x,y,z,atomA,atomB,atomB,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}

void MCS2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	int *note;
	double *x,*y,*z;
	char *Line_End;
	char Line_Date[100];
	FILE *fp;

	if(A > N)
	{
		printf("原子个数超过总原子数量");
		return;
	}

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	Line_End = StoragePath(shape,N,A,N-A,atomA,atomB,atomB,Output);
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	//判断是否存在result.txt，存在则继续实验，不存在则到Cood文件夹读取坐标然后产生随机排布。
	if((fp=fopen(Line_Date,"r"))==NULL){
		printf("\nfile   not   exist\n");
		ReadCood(shape,N,x,y,z);
		srand((unsigned)time(NULL));
		L0NoteInt2(note,N,A,x,y,z);
	} else     
	{
		printf("\nfile   exist\n");  
		ReadFile(Line_Date,note,x,y,z,N);
		fclose(fp);
	}
	
	MCS3_Start(shape,N,A,N-A,note,x,y,z,atomA,atomB,atomB,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}



//蒙特卡洛算法 三合金
void MCS3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int *note;
	double *x,*y,*z;
	char *Line_End;
	char Line_Date[100];
	FILE *fp;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	//判断是否存在result.txt，存在则继续实验，不存在则到Cood文件夹读取坐标然后产生随机排布。
	if((fp=fopen(Line_Date,"r"))==NULL){
		printf("\nfile   not   exist\n");
		ReadCood(shape,N,x,y,z);
		srand((unsigned)time(NULL));
		MixNoteInt3(note,N,A,B);
	//	ShellByShellNoteInt2(note,N,x,y,z);
	} else     
	{
		printf("\nfile   exist\n");  
		ReadFile(Line_Date,note,x,y,z,N);
		fclose(fp);
	}
	
	MCS3_Start(shape,N,A,B,note,x,y,z,atomA,atomB,atomC,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
//	freeSetEnergyNewModel();
}

void MCS3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int i;
	int A;
	int *note;
	int *shell;
	double *x,*y,*z;
	char *Line_End;
	char Line_Date[100];
	FILE *fp;

	shell = Shell_Shape(shape,N);

	A =0;

	for(i = 0;i < N;i ++)
	{
		if(shell[i] > surfaceLayer)
			A ++;
	}
	if(B > N - A)
	{
		printf("B原子超过最大可用原子数\n");
		printf("Core:%d\tSurface:%d\n",A,N -A);
		return;
	}

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	//判断是否存在result.txt，存在则继续实验，不存在则到Cood文件夹读取坐标然后产生随机排布。
	if((fp=fopen(Line_Date,"r"))==NULL){
		printf("\nfile   not   exist\n");
		ReadCood(shape,N,x,y,z);
		srand((unsigned)time(NULL));
		CoreSurfaceNote3(note,N,x,y,z,surfaceLayer,B);
	} else     
	{
		printf("\nfile   exist\n");  
		ReadFile(Line_Date,note,x,y,z,N);
		fclose(fp);
	}
	
	MCS3_Start(shape,N,A,B,note,x,y,z,atomA,atomB,atomC,energy,RSTEP,Output);

	free(shell);
	free(note);
	free(x);
	free(y);
	free(z);
}

void MCS3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int *note;
	double *x,*y,*z;
	char *Line_End;
	char Line_Date[100];
	FILE *fp;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	//判断是否存在result.txt，存在则继续实验，不存在则到Cood文件夹读取坐标然后初始产生。
	if((fp=fopen(Line_Date,"r"))==NULL){
		printf("\nfile   not   exist\n");
		ReadCood(shape,N,x,y,z);
		srand((unsigned)time(NULL));
		FromCoreNoteInt3(note,N,A,B,x,y,z);
		ReadCood(shape,N,x,y,z);
	} else     
	{
		printf("\nfile   exist\n");  
		ReadFile(Line_Date,note,x,y,z,N);
		fclose(fp);
	}
	
	MCS3_Start(shape,N,A,B,note,x,y,z,atomA,atomB,atomC,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}



void MCS3_Start(char *shape,int N,int A,int B,int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int i;
	double a0;
	int clocks = 0;
	int IJKL = 0;
	double *R;
	int numberA,numberB,numberC;
	double biliA,biliB,biliC;
	char *nameA,*nameB,*nameC;
	static char *Line_End;
	char Line_Date[200];
	double E0,E1;
	FILE *fp;
	int tig;
	int MMJ,NNJ;
	int MJ,NJ;

	srand((unsigned)time(NULL));			// 随机种子
	numberA = 0;numberB =0;numberC = 0;
	for(i=0;i<N;i++){
		if(note[i]==0)
			numberA ++;
		else if(note[i]==1)
			numberB ++;
		else
			numberC ++;
	}
	if(numberA != A || numberB != B)
	{
		printf("原子数量不匹配");
		return;
	}
	
	a0 = getLatticeParameter(atomA,atomB,atomC);
	biliA = (double)numberA / N;
	biliB = (double)numberB / N;
	biliC = (double)numberC / N;

	nameA = GetAtomPara(atomA).name; nameB = GetAtomPara(atomB).name;nameC = GetAtomPara(atomC).name;
	R = calloc(N*N,sizeof(double));
	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
	

	biliA = (double)A/N;
	biliB = (double)B/N; 

	Distance(x,y,z,R,N);
	//晶格常数设置为第一种原子
	for(i=0;i<N;i++){
		x[i] = x[i]*a0/(Rmin*sqrt(2));
		y[i] = y[i]*a0/(Rmin*sqrt(2));
		z[i] = z[i]*a0/(Rmin*sqrt(2));
	} 
	Distance(x,y,z,R,N);

		//计算势能数据
	SetEnergyPow(atomA,atomB,atomC);
	setupJohnson();
//	setEnergyNewModel(atomA,atomB,R,N);
		
//	E0 = QSCEnergy3(note,R,atomA,atomB,atomC,N,a0);
//	E0 = QSCCutEnergy3(note,R,atomA,atomB,atomC,N,a0);
//	E0 = JohnsonEnergyCut(note,R,atomA,atomB,atomC,N,a0);
//	E0 = AEAMEnergy(note,R,N);
//	E0 = QSCEnergyNewModel(note,atomA,atomB,N);
	E0 = GetCutEnergyFunction(energy)(note,R,atomA,atomB,atomC,N,a0);

	printf("\nInit MC...\n");
	printf("Shape=%s\tN=%d\n",shape,N);
	printf("E0=%lf\n",E0);
	printf("biliA=%.3f\tbiliB=%.3f\tbiliC=%.3f\n",biliA,biliB,biliC);
	printf("%s=%d\t%s=%d\t%s=%d\n",nameA,numberA,nameB,numberB,nameC,numberC);
	printf("a0=%lf\n",Rmin*sqrt(2));

	E1 = E0;
	printf("\nStart MC...\n");
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\energy.txt");
	fp = fopen(Line_Date,"w");
	fprintf(fp,"\t%d\t%.3lf\n",N,biliA);
	fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);
	while(1)
	{
		tig = 0;
		clocks++;
		if(clocks > RSTEP)
			break;
		MMJ = (int)(N*rand()/(RAND_MAX+1.0));
		NNJ = (int)(N*rand()/(RAND_MAX+1.0));
		MJ = note[MMJ];
		NJ = note[NNJ];
		while(MJ==NJ)                               // 确保随机产生的坐标位置上的原子不同 
		{
			MMJ = (int)(N*rand()/(RAND_MAX+1.0));
			NNJ = (int)(N*rand()/(RAND_MAX+1.0));
			MJ = note[MMJ];
			NJ = note[NNJ];
		}

		note[MMJ] = NJ;
		note[NNJ] = MJ;

	//	E1 = QSCEnergy3(note,R,atomA,atomB,atomC,N,a0);
	//	E1 = QSCCutEnergy3(note,R,atomA,atomB,atomC,N,a0);
	//	E1 = JohnsonEnergyCut(note,R,atomA,atomB,atomC,N,a0);
	//	E1 = AEAMEnergy(note,R,N);
	//	E1 = QSCEnergyNewModel(note,atomA,atomB,N);
		E1 = GetCutEnergyFunction(energy)(note,R,atomA,atomB,atomC,N,a0);

		if(E1 <= E0)								 // 交换后能量更低则保留，更改则变回原值 
		{
			IJKL++;
			tig = 1;
			E0 = E1;
		}
		else
		{
			note[MMJ] = MJ;
			note[NNJ] = NJ;
			tig = 0;
		}

		if(clocks%1000 ==0 )
		{
			printf("clocks=%d\tIJKL=%d\n",clocks,IJKL);
			printf("  E   =%f\n",E0);
		}

		if(tig == 1)
		{
			fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);
			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\result.txt");			
			printResult3(note,x,y,z,N,0,Line_Date);

			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\Diamond.txt");
			printDiamond3(note,x,y,z,atomA,atomB,atomC,N,Line_Date);
		}
	}
	fclose(fp);
	free(R);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");			
	printResult3(note,x,y,z,N,0,Line_Date);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond3(note,x,y,z,atomA,atomB,atomC,N,Line_Date);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	printData(Line_Date,Line_End,N);

}