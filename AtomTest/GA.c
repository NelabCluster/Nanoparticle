#include "GA.h"

void GAPara_Init(GAPARA *para)
{
	memcpy(para, &_defGAPara, sizeof(*para));
}

void GA2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,GAPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,END);
	AtomNum_Init(&atomNum,A,N-A,END);
	GA_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);
}

void GA3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,GAPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,atomC,END);
	AtomNum_Init(&atomNum,A,B,N-A-B,END);
	GA_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);	
}

void GA_InitWithMixing(char *shape,int N, ATOMNUM *atomNum,ALLOY *alloy, PE energy, GAPARA *para,char *Output)
{
	int i;
	GAInstance instance;
	
	//初始化遗传算法实例
	instance.shape = shape;
	instance.N = N;
	instance.energyType = energy;
	Alloy_Copy(&instance.alloy,alloy);
	AtomNum_Copy(&instance.atomNum,atomNum);
	Cood_Init(&instance.cood,N);
	instance.para = *para;
	instance.pop = calloc(para->popSize,sizeof(INDIVIDUAL));
	for( i = 0; i < para->popSize; i++)
		GAIndividual_Init(&instance.pop[i],N);
	GAIndividual_Init(&instance.best,N);
	GAIndividual_Init(&instance.worst,N);
	instance.clocks = 0;
	instance.IJKL = 0;

	//根据构型和原子数读取坐标
	ReadCood1(shape,N,&instance.cood);
	
	//如果需要排序坐标则从里到外排序坐标
	if(instance.para.needOrderCood == 1)
		orderCoodFromCore(instance.cood.x,instance.cood.y,instance.cood.z,instance.N);
	
	//设置随机种子，不然每次启动随机值都一样，每个个体随机产生初始解
	srand((unsigned)time(NULL));
	for( i = 0; i < instance.para.popSize;i++)
		MixNoteInt(instance.pop[i].chrom, instance.N, &instance.atomNum);

	//开始进行遗传算法
	GA_Start(&instance,Output);
	
	//释放遗传算法实例中开辟的空间
	GAIndividual_Free(&instance.worst);
	GAIndividual_Free(&instance.best);
	for(i = 0;i < instance.para.popSize;i++)
		GAIndividual_Free(&instance.pop[i]);
	free(instance.pop);
	Cood_Free(&instance.cood);
	AtomNum_Free(&instance.atomNum);
	Alloy_Free(&instance.alloy);
}

void GA_Start(GAInstance *instance, char *output)
{
	FILE *fp = NULL;
	double a0,E0,E1;
	int i, tig = 0, step = 0;
	char *Line_End;

	srand((unsigned)time(NULL));
	instance->dis.R = calloc(instance->N*instance->N,sizeof(double));
	a0 = getLatticeParameter(&instance->alloy);

	Line_End = StoragePath( instance->shape, instance->N, &instance->alloy, &instance->atomNum, output);
	
	Energy_Init( instance->energyType, &instance->alloy );

	Distance1(&instance->cood,&instance->dis);
	for(i=0;i<instance->N;i++){
		instance->cood.x[i] = instance->cood.x[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.y[i] = instance->cood.y[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.z[i] = instance->cood.z[i]*a0/(instance->dis.Rmin*sqrt(2));
	}
	Distance1(&instance->cood,&instance->dis);
	
	// 计算种群中每个个体的能量
	for( i = 0; i < instance->para.popSize; i++)
		instance->pop[i].energy = GetCutEnergyFunction(instance->energyType)(instance->pop[i].chrom,instance->dis.R,&instance->alloy,instance->N);
	
	// 更新 instance->best 和 instance->worst 
	GA_FindBestAndWorst(instance);
	
	// E0:当前代的最优能量，E1:下一代的最优能量
	// 迭代前 E1 = E0，迭代后 E1 将被更新
	E0 = instance->best.energy;
	E1 = E0;
	
	// 打印遗传算法实例的基本信息
	printf("\nInit GA...\n");
	GA_PrintMsg( instance, Line_End );

	printf("\nStart GA...\n");
	GA_EnergyFile(instance,&fp,Line_End);
	while(1){
		E0 = instance->best.energy;
		instance->clocks++;
	
		//选择
		select_operator(instance);
		//交叉
		crossover_operator(instance);
		//变异
		mutation_operator(instance);
		
		for(i=0;i<instance->para.popSize;i++)
			instance->pop[i].energy = GetCutEnergyFunction(instance->energyType)(instance->pop[i].chrom,instance->dis.R,&instance->alloy,instance->N);
				
		if(instance->para.needAdjustment == 1)
			adjustment_operator(instance);	
	
		//保留最优个体
		changeIndividual(&instance->pop[0],&instance->best,instance->N);

		GA_FindBestAndWorst(instance);
		
		E1 = instance->best.energy;
		if( E1-E0<0.001 && E1-E0>-0.001){
			tig = 0;
			step ++;
		}	else {
			tig = 1;
			step = 0;
			instance->IJKL++;
		}
		if(tig == 1 || (step!=0 && step%100 == 0))
		{	
			E0 = E1;
			printf("clocks=%d\tIJKL=%d\tstep=%d\n",instance->clocks,instance->IJKL,step);
			printf("  E   =%f\n",E0);

			GA_EnergyFile(instance,&fp,Line_End);
			GA_ResultFile(instance,Line_End);	
		}
		

		if( step>instance->para.convergenceGenerations || instance->clocks >= instance->para.maxGenerations )
			break;

	}

	fclose(fp);
	free(instance->dis.R);
	Energy_Free( instance->energyType );

	printf("\nEnd GA...\n");
	GA_PrintMsg( instance, Line_End );
	GA_ResultFile(instance,Line_End);

	printData( Line_End, instance->N );
	
	free(Line_End);
}

void select_operator(GAInstance *instance){
	int i,index,N,POPSIZE;
	double bestEnergy,worstEnergy;
	double p,sum=0;
	double *cvalue;
	struct individual *new_pop,*pop;
	double rou = 3;
	double a1=0;
	double a2 = 1;

	pop = instance->pop;
	N = instance->N;
	POPSIZE = instance->para.popSize;
	bestEnergy = instance->best.energy;
	worstEnergy = instance->worst.energy;

	cvalue = calloc(POPSIZE,sizeof(double));
	new_pop = calloc(POPSIZE,sizeof(struct individual));
	for(i=0;i<POPSIZE;i++)
	{
		pop[i].value = exp(-rou*(pop[i].energy-bestEnergy)/(worstEnergy-bestEnergy+0.0001));
	}

	for(i=0;i<POPSIZE;i++)
		sum+=pop[i].value;

	for(i=0;i<POPSIZE;i++)
		cvalue[i] = pop[i].value/sum;

	for(i=1;i<POPSIZE;i++)
		cvalue[i] = cvalue[i-1]+cvalue[i];

	for(i=0;i<POPSIZE;i++)
	{
		new_pop[i].chrom = calloc(N,sizeof(int));
		p = rand()%1000/1000.0;
		index = 0;
		while(p>cvalue[index])
			index++;
		changeIndividual(&new_pop[i],&pop[index],N);
	}

	for(i=0;i<POPSIZE;i++)
		changeIndividual(&pop[i],&new_pop[i],N);

	free(cvalue);
	for(i=0;i<POPSIZE;i++)
		free(new_pop[i].chrom);
	free(new_pop);
}

void crossover_operator(GAInstance *instance){
	int i,j;
	int length;
	int *index;
	int dNum[3];
	int point,point1,point2,point_left,point_right,temp;
	double p;
	int ch;
	struct individual *pop;
	int N,POPSIZE;
	double pc;
	
	pop = instance->pop;
	N = instance->N;
	POPSIZE = instance->para.popSize;
	pc = instance->para.pc;

	index = calloc(POPSIZE,sizeof(int));
	for(i=0;i<POPSIZE;i++)
		index[i] = i;
	for(i=0;i<POPSIZE;i++)
	{
		point = rand()%(POPSIZE-i);
		temp = index[i];
		index[i] = index[point+i];
		index[point+i] = temp;
	}

	for(i=0;i<POPSIZE-1;i+=2)
	{
		p = rand()%1000/1000.0;
		
		if(p<pc)
		{

			dNum[0] = 0; dNum[1] = 0; dNum[2] = 0;
			point_left = (int)(N*rand()/(RAND_MAX+1.0));
			point_right = (int)(N*rand()/(RAND_MAX+1.0));
			
			if(point_right<point_left)
			{
				point = point_left;
				point_left = point_right;
				point_right = point;
			}

			for(j=point_left;j<=point_right;j++)
			{
				(pop[index[i]].chrom[j] == 0)?(dNum[0]--):(dNum[0]);
				(pop[index[i]].chrom[j] == 1)?(dNum[1]--):(dNum[1]);
				(pop[index[i]].chrom[j] == 2)?(dNum[2]--):(dNum[2]);
				(pop[index[i+1]].chrom[j] == 0)?(dNum[0]++):(dNum[0]);
				(pop[index[i+1]].chrom[j] == 1)?(dNum[1]++):(dNum[1]);
				(pop[index[i+1]].chrom[j] == 2)?(dNum[2]++):(dNum[2]);
						
				ch = pop[index[i]].chrom[j];
				pop[index[i]].chrom[j] = pop[index[i+1]].chrom[j];
				pop[index[i+1]].chrom[j] = ch;
			}
			
			length = N - (point_right - point_left +1);
			for(j=0;j<3;j++)
			{	
				while(dNum[j]>0)
				{
					point = rand()%length;
					if(point>=point_left)
						point1 = point_right+(point+1-point_left);
					else
						point1 = point;

					while(pop[index[i]].chrom[point1] != j)
					{
						point = (point+1)%length;
						if(point>=point_left)
							point1 = point_right+(point+1-point_left);
						else
							point1 = point;
					}

					point = rand()%length;
					if(point>=point_left)
						point2 = point_right+(point+1-point_left);
					else
						point2 = point;

					temp = pop[index[i+1]].chrom[point2];
					while(dNum[temp]>=0)
					{
						point = (point+1)%length;
						if(point>=point_left)
							point2 = point_right+(point+1-point_left);
						else
							point2 = point;

						temp = pop[index[i+1]].chrom[point2]; 
					}

					pop[index[i]].chrom[point1] = temp;
					pop[index[i+1]].chrom[point2] = j;

					dNum[j]--;
					dNum[temp]++;
				}
			}
		}
	}

	free(index);
}

void mutation_operator(GAInstance *instance){
	int i,j,point,ch,N,POPSIZE;
	double p,pm;
	struct individual *pop;

	pop = instance->pop;
	N = instance->N;
	POPSIZE = instance->para.popSize;
	pm = instance->para.pm;

	for(i=0;i<POPSIZE;i++)
		for(j=0;j<N;j++)
		{
			p = rand()%1000/1000.0;
			if(p<pm)
			{
				ch = pop[i].chrom[j];
				point = (int)(N*rand()/(RAND_MAX+1.0));
				while(pop[i].chrom[point] == ch)
				{
					point++;
					if(point==N)
					point = 0;
				}	
				pop[i].chrom[j] = pop[i].chrom[point];
				pop[i].chrom[point] = ch;	
			}
		}
}

void adjustment_operator(GAInstance *instance)
{
	int i,randN1,randN2,tempChrom,POPSIZE,N;
	double tempE,r,rate;
	struct individual *pbest,*pop;
	
	POPSIZE = instance->para.popSize;
	rate = instance->para.adjustmentRate;
	pop = instance->pop;
	N = instance->N;


	for(i=0;i<POPSIZE;i++)
	{
		r = RAND1;
		if(r < rate)
		{
			pbest = pop + i;
			randN1 = (int)(N * RAND1);
			randN2 = (int)(N * RAND1);
			while(pbest->chrom[randN1] == pbest->chrom[randN2])
			{
				randN1 = (int)(N * RAND1);
				randN2 = (int)(N * RAND1);
			}
			tempChrom = pbest->chrom[randN1];
			pbest->chrom[randN1] = pbest->chrom[randN2];
			pbest->chrom[randN2] = tempChrom;
			tempE = GetCutEnergyFunction(instance->energyType)(pbest->chrom,instance->dis.R,&instance->alloy,N); 
			if(tempE <= pbest->energy)
			{
				pbest->energy = tempE;
				continue;
			} 
			tempChrom = pbest->chrom[randN1];
			pbest->chrom[randN1] = pbest->chrom[randN2];
			pbest->chrom[randN2] = tempChrom;
		}
	}
}

void GA_FindBestAndWorst(GAInstance *instance)
{
	int i;
	changeIndividual(&instance->best,&(instance->pop[0]),instance->N);
	changeIndividual(&instance->worst,&(instance->pop[0]),instance->N);
	for(i=0;i<instance->para.popSize;i++)
	{
		if(instance->pop[i].energy<instance->best.energy)
			changeIndividual(&instance->best,&(instance->pop[i]),instance->N);
		if(instance->pop[i].energy>instance->worst.energy)
			changeIndividual(&instance->worst,&(instance->pop[i]),instance->N);
	}
}

void GA_PrintMsg( GAInstance *instance,  char *output )
{	
	int i,tempNum;
	ATOM tempATOM;
	FILE *fp;
	char Line_Date[200];
	
	strcpy( Line_Date, output );
	strcat( Line_Date, "\\GA.txt" );

	fp = fopen( Line_Date, "w" );
	printf( "shape=%s\tN=%d\n", instance->shape, instance->N );
	printf("POPSIZE=%d\tpc=%lf\tpm=%lf\n",instance->para.popSize,instance->para.pc,instance->para.pm);
	fprintf( fp, "shape=%s\tN=%d\n", instance->shape, instance->N );
	fprintf( fp, "POPSIZE=%d\tpc=%lf\tpm=%lf\n",instance->para.popSize,instance->para.pc,instance->para.pm );

	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%.3f\t",GetAtomPara(tempATOM).name,(double)tempNum / instance->N);
		fprintf( fp, "%s=%.3f\t", GetAtomPara(tempATOM).name, (double)tempNum / instance->N );
	}
	printf( "\n" );
	fprintf( fp, "\n" );

	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%d\t",GetAtomPara(tempATOM).name,tempNum);
		fprintf( fp, "%s=%d\t",GetAtomPara(tempATOM).name,tempNum );
	}
	printf("\n");
	printf("best=%lf\tworst=%lf\n",instance->best.energy,instance->worst.energy);
	fprintf( fp, "\n" );
	fprintf( fp, "best=%lf\tworst=%lf\n",instance->best.energy,instance->worst.energy );

	fclose( fp );
}

void GA_EnergyFile(GAInstance *instance, FILE **fp, char *output)
{
	char Line_Date[200];
	
	if ( *fp == NULL )
	{
		strcpy(Line_Date,output);
		strcat(Line_Date,"\\energy.txt");
		*fp = fopen(Line_Date,"w");
		fprintf(*fp,"\t%d\n",instance->N);
	}
	fprintf(*fp,"%6d\t%6d\t%lf\n",instance->clocks,instance->IJKL,instance->best.energy);
}

void GA_ResultFile(GAInstance *instance, char *output)
{
	char Line_Date[200];
	
	//生成结果文件，文件路径名：[output]\\result.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\result.txt");			
	printResult(instance->best.chrom, instance->N, &instance->cood, Line_Date);
	
	//生成绘图文件，文件路径名：[output]\\Diamond.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond(instance->best.chrom, instance->N, &instance->cood, &instance->alloy, Line_Date);
}

void GAIndividual_Init(INDIVIDUAL *one,int N) 
{	
	one->chrom = calloc(N,sizeof(int)); 
}

void changeIndividual(struct individual *to,struct individual *from,int N){
	int i;
	to->energy = from->energy;
	to->number = from->number;
	to->value = from->value;
	for(i=0;i<N;i++){
		to->chrom[i] = from->chrom[i];
	}
}

void GAIndividual_Free(INDIVIDUAL *one) 
{
	one->energy =0; one->value = 0; 
	free(one->chrom); one->chrom = NULL; 
}

void GA3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output)
{
	int i;
	int A;
	int *shell;
	double *x,*y,*z;
	struct individual *pop;

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

	ReadCood(shape,N,x,y,z);
	orderCoodFromCore(x,y,z,N);

	pop = calloc(POPSIZE,sizeof(INDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<POPSIZE;i++){
		pop[i].number = i;
		pop[i].chrom = calloc(N,sizeof(int));
		CoreSurfaceNote3(pop[i].chrom,N,x,y,z,surfaceLayer,B);
	}

	//	GA3_Start(shape,N,A,B,pop,x,y,z,atomA,atomB,atomC,energy,POPSIZE,pm,pc,rate,Output);

	free(shell);
	free(x);
	free(y);
	free(z);
	for(i = 0;i < POPSIZE;i++)
		free(pop[i].chrom);
}

void GA3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output)
{
	int i;
	double *x,*y,*z;
	struct individual *pop;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));

	ReadCood(shape,N,x,y,z);
	orderCoodFromCore(x,y,z,N);

	pop = calloc(POPSIZE,sizeof(INDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<POPSIZE;i++){
		pop[i].number = i;
		pop[i].chrom = calloc(N,sizeof(int));
		FromCoreNoteInt3(pop[i].chrom,N,A,B,x,y,z);
	}

	//	GA3_Start(shape,N,A,B,pop,x,y,z,atomA,atomB,atomC,energy,POPSIZE,pm,pc,rate,Output);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < POPSIZE;i++)
		free(pop[i].chrom);
}

