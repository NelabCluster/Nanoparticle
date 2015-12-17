#include "PSO.h"

void PSO2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output,PSOIMPROVE improve)
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
		MixNoteInt3(pop[i].chrom,N,A,N-A);
		pop[i].swapList = initSwapSequence(10,N);
		pop[i].pbest = calloc(1,sizeof(PSOINDIVIDUAL));
		pop[i].pbest->chrom = calloc(N,sizeof(int));
	}

	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,Output,improve);

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

void PSO2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output)
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

	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,Output,PSOYesImprove);

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


void PSO2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output)
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

	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,Output,PSOYesImprove);

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

void PSO2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *Output)
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

	PSO2_Start(shape,N,A,pop,x,y,z,atomA,atomB,energy,Popsize,w,rate,Output,PSOYesImprove);

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

void PSO2_Start(char *shape,int N,int A,PSOINDIVIDUAL *pop,double *x,double *y,double *z,ATOM atomA,ATOM atomB,PE energy,int POPSIZE,double w,double rate,char *Output,PSOIMPROVE improve)
{
	int i;
	FILE *fp;
	double E0,E1;
	int clocks = 0;
	int IJKL = 0;
	int tig = 0;
	int step = 0;
	char *Line_End;
	char Line_Date[100];
	double *R;
	PSOINDIVIDUAL *best_pop;
	double a0;
	int B;
	double biliA,biliB;
	char *nameA,*nameB;

	srand((unsigned)time(NULL));
	B = N - A;
	R = calloc(N*N,sizeof(double));
	a0 = getLatticeParameter(atomA,atomB,atomB);
	nameA = GetAtomPara(atomA).name; nameB = GetAtomPara(atomB).name;
	biliA = (double)A / N;
	biliB = (double)B / N;
	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomB,Output);

	Distance(x,y,z,R,N);
	for(i=0;i<N;i++){
		x[i] = x[i]*a0/(Rmin*sqrt(2));
		y[i] = y[i]*a0/(Rmin*sqrt(2));
		z[i] = z[i]*a0/(Rmin*sqrt(2));
	}
	Distance(x,y,z,R,N);

	SetEnergyPow(atomA,atomB,atomB);
	setupJohnson();
//	setEnergyNewModel(atomA,atomB,R,N);

	//初始化
	for(i=0;i<POPSIZE;i++){
		pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomB,N,a0);
	//	pop[i].energy = QSCEnergyNewModel(pop[i].chrom,atomA,atomB,N);
		changePSOIndividual(pop[i].pbest,&(pop[i]),N);
	}
	
	best_pop = pop[0].pbest;
	updateGbest(&best_pop,pop,POPSIZE);
	
	E0 = best_pop->energy;
	E1 = E0;
	
	printf("\nInit PSO...\n");
	printf("shape=%s\tN=%d\n",shape,N);
	printf("POPSIZE=%d\tw=%lf\trate=%lf\n",POPSIZE,w,rate);
	printf("biliA=%.3f\tbiliB=%.3f\n",biliA,biliB);
	printf("%s=%d\t%s=%d\t\n",nameA,A,nameB,B);	
	printf("best=%lf\n",best_pop->energy);


	printf("\nStart PSO...\n");
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\energy.txt");
	fp = fopen(Line_Date,"w");
	fprintf(fp,"\t%d\t%3.1lf\n",N,biliA);
	fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);

	while(1)
	{
		clocks++;

		//更新交换序
		updateSwapList(w,pop,best_pop,POPSIZE,N);

		//更新位置
		upDateChorm(pop, POPSIZE);

		//更新pbest
		updatePbest(R, pop,atomA,atomB,energy,POPSIZE,N);
		
		//位置调整
		if(improve == PSOYesImprove)
			PSOMutation(rate,R,pop,atomA,atomB,energy,POPSIZE,N);

		//更新gbest
		updateGbest(&best_pop,pop,POPSIZE);

		E1 = best_pop->energy;
		if( E1-E0<0.001 && E1-E0>-0.001){
			tig = 0;
			step ++;
		}	else {
			tig = 1;
			step = 0;
			IJKL++;
		}
		if(tig == 1 || (step!=0 && step % 100 == 0))
		{	
			E0 = E1;
			printf("clocks=%d\tIJKL=%d\tstep=%d\n",clocks,IJKL,step);
			printf("  E   =%f\n",E0);
			fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);

			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\result.txt");			
			printResult3(best_pop->chrom,x,y,z,N,0,Line_Date);

			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\Diamond.txt");
			printDiamond3(best_pop->chrom,x,y,z,atomA,atomB,atomB,N,Line_Date);
		}

		if(step>1000)
			break;

	}
	fclose(fp);
	free(R);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");			
	printResult3(best_pop->chrom,x,y,z,N,0,Line_Date);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond3(best_pop->chrom,x,y,z,atomA,atomB,atomB,N,Line_Date);
	
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	printData(Line_Date,Line_End,N);
	
//	freeSetEnergyNewModel();
}

void updateSwapList(double w, PSOINDIVIDUAL *pop, PSOINDIVIDUAL *gbest, int POPSIZE, int N)
{
	int i;
	PSOSWAP *swapList1,*swapList2;
	for(i = 0; i < POPSIZE; i++)
	{
		swapList1 = swapSubtraction(gbest,&pop[i],N);
		swapList2 = swapSubtraction(pop[i].pbest,&pop[i],N);
		pop[i].swapList = swapAddition(w,&pop[i],swapList1,swapList2);
	}
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

void upDateChorm(PSOINDIVIDUAL *pop, int POPSIZE)
{
	int i;
	int temp0,temp1;
	PSOINDIVIDUAL *p;
	PSOSWAP *s;
	
	for(i=0;i<POPSIZE;i++)
	{
		p = &pop[i];
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

void PSOMutation(double rate, double *R, PSOINDIVIDUAL *pop, ATOM atomA, ATOM atomB, PE energy, int POPSIZE, int N)
{
	int i,randN1,randN2,tempChrom;
	double tempE,r;
	PSOINDIVIDUAL *pbest;

	for(i=0;i<POPSIZE;i++)
	{
		r = RAND1;
		if(r < rate)
		{
			pbest = pop[i].pbest;
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
			tempE = GetCutEnergyFunction(energy)(pbest->chrom,R,atomA,atomB,atomB,N,getLatticeParameter(atomA,atomB,atomB)); 
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

void updatePbest(double *R, PSOINDIVIDUAL *pop, ATOM atomA, ATOM atomB, PE energy, int POPSIZE, int N)
{
	int i;

	for(i=0;i<POPSIZE;i++)
	{
		pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomB,N,getLatticeParameter(atomA,atomB,atomB));
		if(pop[i].energy <= pop[i].pbest->energy)
		{
			changePSOIndividual(pop[i].pbest,&pop[i],N);
		}
	}
}

void updateGbest(PSOINDIVIDUAL **best_pop,PSOINDIVIDUAL *pop,int POPSIZE)
{
	int i;
	
	for(i =0 ;i< POPSIZE;i++)
	{
		if(pop[i].pbest->energy <= (*best_pop)->energy)
		{
			*best_pop = pop[i].pbest;
		}
	}
//	changePSOIndividual(best_pop,pop[index].pbest,N);
}



void changePSOIndividual(PSOINDIVIDUAL *one,PSOINDIVIDUAL *two,int N){
	int i;
//	PSOSWAP *swap1,*swap2;
	one->energy=two->energy;

	for(i=0;i<N;i++){
		one->chrom[i] = two->chrom[i];
	}
	
/*	swap1 = one->swapList;
	while(swap1 != NULL){
		swap2 = swap1;
		swap1 = swap1->next;
		free(swap2);
	}
	one->swapList = NULL;

	swap1 = one->swapList;
	swap2 = two->swapList;
	while(swap2 != NULL){
		swap1 = creatSwap(swap2->element[0],swap2->element[1]);
		swap1 = swap1->next;
		swap2 = swap2->next;
	}*/
}
