#include "PSO.h"

void PSOPara_Init(PSOPARA *para)
{
	memcpy(para, &_defPSOPara, sizeof(*para));
}

void PSO2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,PSOPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,END);
	AtomNum_Init(&atomNum,A,N-A,END);
	PSO_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);
}

void PSO3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,PSOPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,atomC,END);
	AtomNum_Init(&atomNum,A,B,N-A-B,END);
	PSO_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);	
}

void PSO_InitWithMixing(char *shape,int N,ATOMNUM *atomNum,ALLOY *alloy,PE energy,PSOPARA *para,char *output)
{
	int i;
	PSOInstance instance;
	
	// 初始化粒子群实例
	instance.shape = shape;
	instance.N = N;
	instance.energyType = energy;
	Alloy_Copy(&instance.alloy,alloy);
	AtomNum_Copy(&instance.atomNum,atomNum);
	Cood_Init(&instance.cood,N);
	instance.para = *para;
	instance.pop = calloc(para->popSize,sizeof(PSOINDIVIDUAL));
	for( i = 0; i < para->popSize; i++)
		PSOIndividual_Init(&instance.pop[i],N);
	instance.gbest = NULL;
	instance.clocks = 0;
	instance.IJKL = 0;

	//根据构型和原子数读取坐标
	ReadCood1(shape,N,&instance.cood);
	
	//如果需要排序坐标则从里到外排序坐标
	if(instance.para.needOrderCood == 1)
		orderCoodFromCore(instance.cood.x,instance.cood.y,instance.cood.z,instance.N);
	
	//设置随机种子，不然每次启动随机值都一样
	//每个个体随机产生初始解 和 10组交换子
	srand((unsigned)time(NULL));
	for( i = 0; i < instance.para.popSize;i++)
	{
		MixNoteInt(instance.pop[i].chrom, instance.N, &instance.atomNum);
		instance.pop[i].swapList = initSwapSequence(10,N);
	}

	//开始进行遗传算法
	PSO_Start(&instance,output);
	
	//释放遗传算法实例中开辟的空间
	for(i = 0;i < instance.para.popSize;i++)
		PSOIndividual_Free(&instance.pop[i]);
	free(instance.pop);
	Cood_Free(&instance.cood);
	AtomNum_Free(&instance.atomNum);
	Alloy_Free(&instance.alloy);
}

void PSO_Start(PSOInstance *instance,char *output)
{
	int i;
	FILE *fp = NULL;
	double E0,E1;
	int tig = 0;
	int step = 0;
	char *Line_End;
	double a0;

	srand((unsigned)time(NULL));
	instance->dis.R = calloc( instance->N*instance->N,sizeof(double));
	a0 = getLatticeParameter( &instance->alloy );

	Line_End = StoragePath( instance->shape, instance->N, &instance->alloy, &instance->atomNum, output );
	
	Energy_Init( instance->energyType, &instance->alloy );

	Distance1(&instance->cood,&instance->dis);
	for(i=0;i<instance->N;i++){
		instance->cood.x[i] = instance->cood.x[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.y[i] = instance->cood.y[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.z[i] = instance->cood.z[i]*a0/(instance->dis.Rmin*sqrt(2));
	}
	Distance1(&instance->cood,&instance->dis);

	//初始化
	for( i = 0; i < instance->para.popSize; i++)
		instance->pop[i].energy = GetCutEnergyFunction1(instance->energyType)(instance->pop[i].chrom,instance->dis.R,&instance->alloy,instance->N);
	updatePbest( instance );
	updateGbest( instance );
	
	E0 = instance->gbest->energy;
	
	printf( "\nInit PSO...\n" );
	PSO_PrintMsg( instance, Line_End );

	printf( "\nStart PSO...\n" );
	PSO_EnergyFile( instance, &fp, Line_End );

	while(1)
	{
		instance->clocks++;

		//更新交换序
		updateSwapList( instance );

		//更新位置
		upDateChorm( instance );

		for( i = 0; i < instance->para.popSize; i++)
			instance->pop[i].energy = GetCutEnergyFunction1(instance->energyType)(instance->pop[i].chrom,instance->dis.R,&instance->alloy,instance->N);

		//更新pbest
		updatePbest( instance );

		//位置调整
		if( instance->para.needAdjustment == 1 )
			PSOMutation( instance );

		//更新gbest
		updateGbest( instance );

		E1 = instance->gbest->energy;
		if( E1-E0<0.001 && E1-E0>-0.001){
			tig = 0;
			step ++;
		}	else {
			tig = 1;
			step = 0;
			instance->IJKL++;
		}
		if(tig == 1 || (step!=0 && step % 100 == 0))
		{	
			E0 = E1;
			printf( "clocks=%d\tIJKL=%d\tstep=%d\n", instance->clocks, instance->IJKL, step );
			printf( "  E   =%f\n", E0 );

			PSO_EnergyFile( instance, &fp, Line_End );
			PSO_ResultFile( instance, Line_End );
		}

		if( step > instance->para.convergenceGenerations || instance->clocks > instance->para.maxGenerations )	break;
	}

	fclose(fp);
	free( instance->dis.R );
	
	Energy_Free( instance->energyType );

	PSO_PrintMsg( instance, Line_End );
	PSO_ResultFile( instance, Line_End );
	
	printData( Line_End, instance->N );

	free(Line_End);
}

void updatePbest( PSOInstance *instance )
{
	int i;

	for( i = 0; i < instance->para.popSize; i++ )
	{
		if( instance->pop[i].energy <= instance->pop[i].pbest->energy)
			changePSOIndividual( instance->pop[i].pbest, &instance->pop[i], instance->N );
	}
}

void updateGbest( PSOInstance *instance )
{
	int i;
	
	if( instance->gbest == NULL ) instance->gbest = instance->pop[0].pbest;
	
	for( i = 0 ; i < instance->para.popSize; i++ )
	{
		if( instance->pop[i].pbest->energy <= instance->gbest->energy)
		{
			instance->gbest = instance->pop[i].pbest;
		}
	}
}

void updateSwapList( PSOInstance *instance )
{
	int i;
	PSOSWAP *swapList1,*swapList2;

	for( i = 0; i < instance->para.popSize; i++ )
	{
		swapList1 = swapSubtraction( instance->gbest, &instance->pop[i], instance->N );
		swapList2 = swapSubtraction( instance->pop[i].pbest, &instance->pop[i], instance->N );
		instance->pop[i].swapList = swapAddition( instance->para.w, &instance->pop[i], swapList1, swapList2 );
	}
}

void upDateChorm( PSOInstance *instance )
{
	int i;
	int temp0,temp1;
	PSOINDIVIDUAL *p;
	PSOSWAP *s;
	
	for( i = 0; i < instance->para.popSize; i++ )
	{
		p = &instance->pop[i];
		s = p->swapList->next;
		while(s != NULL)
		{
			temp0 = p->chrom[s->element[0]];
			temp1 = p->chrom[s->element[1]];
			if(temp0 != temp1)
			{
				p->chrom[s->element[0]] = temp1;
				p->chrom[s->element[1]] = temp0;
			}
			s = s->next;
		}
	}
}

void PSOMutation( PSOInstance *instance )
{
	int i,randN1,randN2,tempChrom;
	double tempE,r;
	PSOINDIVIDUAL *pbest;

	for( i = 0; i < instance->para.popSize; i++ )
	{
		r = RAND1;
		if( r < instance->para.rate )
		{
			pbest = instance->pop[i].pbest;
			randN1 = RANDINT( instance->N );
			randN2 = RANDINT( instance->N );
			while( pbest->chrom[randN1] == pbest->chrom[randN2])
			{
				randN1 = RANDINT( instance->N );
				randN2 = RANDINT( instance->N );
			}
			tempChrom = pbest->chrom[randN1];
			pbest->chrom[randN1] = pbest->chrom[randN2];
			pbest->chrom[randN2] = tempChrom;
			tempE = GetCutEnergyFunction1(instance->energyType)(pbest->chrom,instance->dis.R,&instance->alloy,instance->N); 
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

void PSO_PrintMsg(PSOInstance *instance, char *output)
{
	int i,tempNum;
	ATOM tempATOM;
	FILE *fp;
	char Line_Date[200];
	
	strcpy( Line_Date, output );
	strcat( Line_Date, "\\PSO.txt" );

	fp = fopen( Line_Date, "w" );
	printf( "shape=%s\tN=%d\n", instance->shape, instance->N );
	printf( "POPSIZE=%d\tw=%lf\trate=%lf\n", instance->para.popSize, instance->para.w, instance->para.rate );
	fprintf( fp, "shape=%s\tN=%d\n", instance->shape, instance->N );
	fprintf( fp, "POPSIZE=%d\tw=%lf\trate=%lf\n", instance->para.popSize, instance->para.w, instance->para.rate );

	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%.3f\t",GetAtomPara(tempATOM).name,(double)tempNum / instance->N);
		fprintf( fp, "%s=%.3f\t",GetAtomPara(tempATOM).name,(double)tempNum / instance->N);
	}
	printf("\n");
	fprintf( fp, "\n");
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf("%s=%d\t",GetAtomPara(tempATOM).name,tempNum);
		fprintf(fp,"%s=%d\t",GetAtomPara(tempATOM).name,tempNum);
	}
	printf("\n");
	printf("best=%lf\n",instance->gbest->energy);
	fprintf(fp,"\n");
	fprintf(fp,"best=%lf\n",instance->gbest->energy);
}

void PSO_EnergyFile( PSOInstance *instance, FILE **fp, char *output )
{
	char Line_Date[200];

	if ( *fp == NULL )
	{
		strcpy(Line_Date,output);
		strcat(Line_Date,"\\energy.txt");
		*fp = fopen(Line_Date,"w");
		fprintf(*fp,"\t%d\n",instance->N);
	}
	fprintf(*fp,"%6d\t%6d\t%lf\n",instance->clocks,instance->IJKL,instance->gbest->energy);
}

void PSO_ResultFile( PSOInstance *instance, char *output )
{
	char Line_Date[200];

	//生成结果文件，文件路径名：[output]\\result.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\result.txt");			
	printResult(instance->gbest->chrom, instance->N, &instance->cood, Line_Date);
	
	//生成绘图文件，文件路径名：[output]\\Diamond.txt
	strcpy(Line_Date,output);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond(instance->gbest->chrom, instance->N, &instance->cood, &instance->alloy, Line_Date);
}

void PSOIndividual_Init(PSOINDIVIDUAL *one,int N)
{
	one->chrom = calloc(N,sizeof(int));
	one->energy = 0;
	one->swapList = NULL;
	one->pbest = calloc(1,sizeof(PSOINDIVIDUAL));
	one->pbest->chrom = calloc(N,sizeof(int));
}

void changePSOIndividual(PSOINDIVIDUAL *one,PSOINDIVIDUAL *two,int N)
{
	int i;

	one->energy=two->energy;
	for(i=0;i<N;i++){
		one->chrom[i] = two->chrom[i];
	}
}

void PSOIndividual_Free(PSOINDIVIDUAL *one)
{
	PSOSWAP *swap,*swapNext;

	free(one->chrom);
	one->chrom = NULL;
	one->energy = 0;
	swapNext = one->swapList;
	while(swapNext != NULL)
	{
		swap = swapNext;
		swapNext = swap->next;
		free(swap);
	}
	one->swapList = NULL;
	free(one->pbest->chrom);
	free(one->pbest);
	one->pbest = NULL;
}

void printfPSOIndividual(PSOINDIVIDUAL *individual,int N)
{
	int i;
	printf("Energy:%f\n",individual->energy);
	for(i=0 ;i < N;i ++)
	{
		printf("%3d",individual->chrom[i]);
	}
	printf("\n");
}

PSOSWAP* initSwapSequence(int size,int N)
{
	int i;
	int rand1=0,rand2=1;
	PSOSWAP *tempSwap,*swap,*startSwap;
	
	startSwap = creatSwap(0,0);
	swap = startSwap;
	for(i=0;i<size;i++)
	{
		do{
			rand1 = (int)(N*RAND1);
			rand2 = (int)(N*RAND1);
		}while(rand1==rand2);
		tempSwap = creatSwap(rand1,rand2);
		swap->next = tempSwap;
		swap = swap->next;
	}
	return startSwap;
}

PSOSWAP* creatSwap(int element1,int element2)
{
	PSOSWAP *swap;
	swap = (PSOSWAP *)malloc(sizeof(PSOSWAP));
	swap->element[0] = element1;
	swap->element[1] = element2;
	swap->next = NULL;
	return swap;
}

PSOSWAP* swapSubtraction(PSOINDIVIDUAL *best, PSOINDIVIDUAL *pop, int N)
{
	int i,j;
	int temp;
	int *first,*second;
	PSOSWAP *head,*s;
	
	first = best->chrom;
	second = malloc(N * sizeof(int));
	head = creatSwap(0,0);
	s = head;
	for(i = 0;i <N; i++)
		second[i] = pop->chrom[i];
	for (i = 0; i < N; i++) {
			if (second[i] != first[i]) {
				for (j = i+1; j < N; j++) {
					if (second[j] != second[i]) {
						if (second[j] != first[j]) {
							if(second[j] == first[i]) {
								temp = second[j];
								second[j] = second[i];
								second[i] = temp;
								s->next	= creatSwap(i,j);
								s = s->next;
								break;
							}
						}
					}
				}
			}	
		}

	free(second);
	return head;
}

PSOSWAP* swapAddition(double w, PSOINDIVIDUAL *pop, PSOSWAP *swapList1, PSOSWAP *swapList2)
{
	double r1,r2;
	PSOSWAP *newSwapList,*lastSwap;
	
	newSwapList = creatSwap(0,0);
	w = 0.5;
	r1 = RAND1;
	r2 = RAND1;
	lastSwap = swapMultiplication(w,pop->swapList,newSwapList);
	lastSwap = swapMultiplication(r1,swapList1,lastSwap);
	lastSwap = swapMultiplication(r2,swapList2,lastSwap);

	return newSwapList;
}

PSOSWAP* swapMultiplication(double rate, PSOSWAP *swapList, PSOSWAP *newSwapList)
{
	double r;
	PSOSWAP *p;

	p = swapList;
	swapList = swapList->next;
	free(p);
	while(swapList != NULL)
	{
		r = RAND1;
		if(r < rate)
		{
			newSwapList->next = swapList;
			newSwapList = newSwapList->next;
			swapList = swapList->next;
		} else
		{
			p = swapList;
			swapList = swapList->next;
			free(p);
		}
	}
	newSwapList->next = NULL;
	return newSwapList;
}

void printfPSOSwap(PSOSWAP *head)
{
	int i = 0;
	PSOSWAP *p;
	p = head;
	while(p!=NULL)
	{
		printf("%5d-%d-%d\t",i,p->element[0],p->element[1]);
		p = p->next;
		i++;
	}
	printf("\n");
}

void PSO2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output)
{
	int i;
	double *x,*y,*z;
	PSOINDIVIDUAL *pop;
	PSOSWAP *swap,*swapNext;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));

	ReadCood(shape,N,x,y,z);
	orderCoodFromCore(x,y,z,N);

	pop = calloc(Popsize,sizeof(PSOINDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<Popsize;i++){
		pop[i].chrom = calloc(N,sizeof(int));
		FromCoreNoteInt3(pop[i].chrom,N,A,N-A,x,y,z);
		pop[i].swapList = initSwapSequence(10,N);
		pop[i].pbest = calloc(1,sizeof(PSOINDIVIDUAL));
		pop[i].pbest->chrom = calloc(N,sizeof(int));
	}

	//	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,output,PSOYesImprove);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < Popsize;i++){
		free(pop[i].chrom);
		free(pop[i].pbest->chrom);
		free(pop[i].pbest);
		swapNext = pop[i].swapList;
		while(swapNext != NULL)
		{
			swap = swapNext;
			swapNext = swap->next;
			free(swap);
		}
	}
	free(pop);
}


void PSO2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output)
{
	int i;
	double *x,*y,*z;
	PSOINDIVIDUAL *pop;
	PSOSWAP *swap,*swapNext;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));

	ReadCood(shape,N,x,y,z);
	orderCoodFromCore(x,y,z,N);

	pop = calloc(Popsize,sizeof(PSOINDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<Popsize;i++){
		pop[i].chrom = calloc(N,sizeof(int));
		PhaseSeparationNoteInt2(pop[i].chrom,N,A,x,y,z);
		pop[i].swapList = initSwapSequence(10,N);
		pop[i].pbest = calloc(1,sizeof(PSOINDIVIDUAL));
		pop[i].pbest->chrom = calloc(N,sizeof(int));
	}

	//	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,output,PSOYesImprove);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < Popsize;i++){
		free(pop[i].chrom);
		free(pop[i].pbest->chrom);
		free(pop[i].pbest);
		swapNext = pop[i].swapList;
		while(swapNext != NULL)
		{
			swap = swapNext;
			swapNext = swap->next;
			free(swap);
		}
	}
	free(pop);
}

void PSO2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output)
{
	int i;
	double *x,*y,*z;
	PSOINDIVIDUAL *pop;
	PSOSWAP *swap,*swapNext;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));

	ReadCood(shape,N,x,y,z);
	orderCoodFromCore(x,y,z,N);

	pop = calloc(Popsize,sizeof(PSOINDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<Popsize;i++){
		pop[i].chrom = calloc(N,sizeof(int));
		L0NoteInt2(pop[i].chrom,N,A,x,y,z);
		pop[i].swapList = initSwapSequence(10,N);
		pop[i].pbest = calloc(1,sizeof(PSOINDIVIDUAL));
		pop[i].pbest->chrom = calloc(N,sizeof(int));
	}

	//	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,output,PSOYesImprove);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < Popsize;i++){
		free(pop[i].chrom);
		free(pop[i].pbest->chrom);
		free(pop[i].pbest);
		swapNext = pop[i].swapList;
		while(swapNext != NULL)
		{
			swap = swapNext;
			swapNext = swap->next;
			free(swap);
		}
	}
	free(pop);
}