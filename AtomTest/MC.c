#include "MC.h"

void MCPara_Init(MCPARA *para)
{
	memcpy(para, &_defMCPara, sizeof(*para));
}

void MC2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,MCPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,END);
	AtomNum_Init(&atomNum,A,N-A,END);
	MC_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);	

}

void MC3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,MCPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,atomC,END);
	AtomNum_Init(&atomNum,A,B,N-A-B,END);
	MC_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);
}

void MC_InitWithMixing(char *shape, int N, ATOMNUM *atomNum, ALLOY *alloy, PE energy, MCPARA *para, char *output)
{
	MCInstance instance;
	
	// 初始化蒙特卡罗法算法实例
	instance.shape = shape;
	instance.N = N;
	instance.energyType = energy;
	Alloy_Copy( &instance.alloy, alloy );
	AtomNum_Copy( &instance.atomNum, atomNum );
	Cood_Init( &instance.cood, N );
	instance.para = *para;
	MCIndividual_Init( &instance.one, N );
	instance.clocks = 0;
	instance.IJKL = 0;
	
	//根据构型和原子数读取坐标
	ReadCood1(shape,N,&instance.cood);
	
	//设置随机种子，不然每次启动随机值都一样，个体随机产生初始解
	srand((unsigned)time(NULL));
	MixNoteInt(instance.one.note, instance.N, &instance.atomNum);
	
	//开始进行蒙特卡罗方法
	MC_Start( &instance, output);
	
	//释放蒙特卡罗算法实例中开辟的空间
	MCIndividual_Free( &instance.one );
	Cood_Free( &instance.cood );
	AtomNum_Free( &instance.atomNum );
	Alloy_Free( &instance.alloy );
}

void MC_Start( MCInstance *instance, char *output ){
	int i;
	double a0;
	static char *Line_End;
	double E1;
	FILE *fp = NULL;
	int tig;
	int MMJ,NNJ;
	int MJ,NJ;

	srand((unsigned)time(NULL));			// 随机种子

	a0 = getLatticeParameter(&instance->alloy);

	instance->dis.R = calloc( instance->N*instance->N,sizeof(double));
	Line_End = StoragePath( instance->shape, instance->N, &instance->alloy, &instance->atomNum, output);
	
	Distance1(&instance->cood,&instance->dis);
	for(i=0;i<instance->N;i++){
		instance->cood.x[i] = instance->cood.x[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.y[i] = instance->cood.y[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.z[i] = instance->cood.z[i]*a0/(instance->dis.Rmin*sqrt(2));
	}
	Distance1(&instance->cood,&instance->dis);

	instance->one.E = GetEnergyFunction1(instance->energyType)(instance->one.note,instance->dis.R,instance->alloy.atoms,instance->alloy.atomTypeCount,instance->N);

	printf("\nInit MC...\n");
	MC_PrintMsg(instance);

	printf("\nStart MC...\n");
	MC_EnergyFile(instance,&fp,Line_End);

	while(1)
	{
		tig = 0;
		instance->clocks++;
		if( instance->clocks > instance->para.RSTEP)
			break;
		MMJ = RANDINT( instance->N );
		NNJ = RANDINT( instance->N );
		MJ = instance->one.note[MMJ];
		NJ = instance->one.note[NNJ];
		while(MJ==NJ)                               // 确保随机产生的坐标位置上的原子不同 
		{
			MMJ = RANDINT( instance->N );
			NNJ = RANDINT( instance->N );
			MJ = instance->one.note[MMJ];
			NJ = instance->one.note[NNJ];
		}

		instance->one.note[MMJ] = NJ;
		instance->one.note[NNJ] = MJ;

		E1 = GetEnergyFunction1(instance->energyType)(instance->one.note,instance->dis.R,instance->alloy.atoms,instance->alloy.atomTypeCount,instance->N);

		if(E1 <= instance->one.E)								 // 交换后能量更低则保留，更改则变回原值 
		{
			instance->IJKL++;
			tig = 1;
			instance->one.E = E1;
		}
		else
		{
			instance->one.note[MMJ] = MJ;
			instance->one.note[NNJ] = NJ;
			tig = 0;
		}

		if( instance->clocks%1000 ==0 )
		{
			printf("clocks=%d\tIJKL=%d\n", instance->clocks, instance->IJKL);
			printf("  E   =%f\n",instance->one.E);
		}

		if(tig == 1)
		{
			MC_EnergyFile( instance, &fp, Line_End );
			MC_ResultFile( instance, Line_End );
		}
	}
	fclose(fp);
	free(instance->dis.R);

	MC_ResultFile( instance, Line_End );

	free(Line_End);
}

void MC_PrintMsg(MCInstance *instance)
{
	int i,tempNum;
	ATOM tempATOM;

	printf("Shape=%s\tN=%d\n",instance->shape,instance->N);
	printf("RSTEP=%ld\tT=%d\n",instance->para.RSTEP,instance->para.temperature);
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%.3f\t",GetAtomPara(tempATOM).name,(double)tempNum / instance->N);
	}
	printf("\n");
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%d\t",GetAtomPara(tempATOM).name,tempNum);
	}
	printf("\n");
	printf("E=%lf\n",instance->one.E);
	printf("a0=%lf\n",instance->dis.Rmin*sqrt(2));
}

void MC_EnergyFile(MCInstance *instance, FILE **fp, char *output)
{
	char Line_Date[200];

	if ( *fp == NULL )
	{
		strcpy(Line_Date,output);
		strcat(Line_Date,"\\energy.txt");
		*fp = fopen(Line_Date,"w");
		fprintf(*fp,"\t%d\n",instance->N);
	}
	fprintf(*fp,"%6d\t%6d\t%lf\n",instance->clocks,instance->IJKL,instance->one.E);
}

void MC_ResultFile(MCInstance *instance, char *output)
{
	char Line_Date[200];

	//生成结果文件，文件路径名：[output]\\result.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\result.txt");			
	printResult(instance->one.note, instance->N, &instance->cood, Line_Date);
	
	//生成绘图文件，文件路径名：[output]\\Diamond.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond(instance->one.note, instance->N, &instance->cood, &instance->alloy, Line_Date);
}

void MCIndividual_Init(MCINDIVIDUAL *one,int N)
{
	one->E = 0;
	one->note = calloc(N,sizeof(int));
}

void MCIndividual_Free(MCINDIVIDUAL *one)
{
	one->E = 0;
	free(one->note);
	one->note = NULL;
}

void MC3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int *note;
	double *x,*y,*z;
	char *Line_End = NULL;
	char Line_Date[100];
	FILE *fp;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	note = calloc(N,sizeof(int));

	//	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
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
	
	//MC3_Start(shape,N,A,B,note,x,y,z,atomA,atomB,atomC,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}

void MC2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	int *note;
	double *x,*y,*z;
	char *Line_End = NULL;
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

	//	Line_End = StoragePath(shape,N,A,N-A,atomA,atomB,atomB,Output);
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
	
	//MC3_Start(shape,N,A,N-A,note,x,y,z,atomA,atomB,atomB,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}

void MC2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	int *note;
	double *x,*y,*z;
	char *Line_End = NULL;
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

	//	Line_End = StoragePath(shape,N,A,N-A,atomA,atomB,atomB,Output);
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
	
	//MC3_Start(shape,N,A,N-A,note,x,y,z,atomA,atomB,atomB,energy,RSTEP,Output);

	free(note);
	free(x);
	free(y);
	free(z);
}

void MC2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output)
{
	if(A > N)
	{
		printf("原子个数超过总原子数量");
		return;
	}
	MC3_InitWithFromCore(shape,N,A,N-A,atomA,atomB,atomB,energy,RSTEP,Output);
}

void MC3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output){
	int i;
	int A;
	int *note;
	int *shell;
	double *x,*y,*z;
	char *Line_End = NULL;
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

	//	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);
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
	
	//MC3_Start(shape,N,A,B,note,x,y,z,atomA,atomB,atomC,energy,RSTEP,Output);

	free(shell);
	free(note);
	free(x);
	free(y);
	free(z);
}