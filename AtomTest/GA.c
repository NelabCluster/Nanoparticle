#include "GA.h"

//遗传算法

void GA2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output)
{
	GA3_InitWithMixing(shape,N,A,N-A,atomA,atomB,atomB,energy,POPSIZE,pc,pm,rate,Output);
}

void GA3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output)
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
		MixNoteInt3(pop[i].chrom,N,A,B);
	}

	GA3_Start(shape,N,A,B,pop,x,y,z,atomA,atomB,atomC,energy,POPSIZE,pm,pc,rate,Output);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < POPSIZE;i++)
		free(pop[i].chrom);
	free(pop);
	
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

	GA3_Start(shape,N,A,B,pop,x,y,z,atomA,atomB,atomC,energy,POPSIZE,pm,pc,rate,Output);

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

	GA3_Start(shape,N,A,B,pop,x,y,z,atomA,atomB,atomC,energy,POPSIZE,pm,pc,rate,Output);

	free(x);
	free(y);
	free(z);
	for(i = 0;i < POPSIZE;i++)
		free(pop[i].chrom);
}



void GA3_Start(char *shape,int N,int A,int B,INDIVIDUAL *pop,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pm,double pc,double rate,char *Output){
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
	struct individual best_pop;
	struct individual worst_pop;
	double a0;
	int C;
	double biliA,biliB,biliC;
	char *nameA,*nameB,*nameC;

	srand((unsigned)time(NULL));
	C = N - A - B;
	R = calloc(N*N,sizeof(double));
	a0 = getLatticeParameter3(atomA,atomB,atomC);
	nameA = GetAtomPara(atomA).name; nameB = GetAtomPara(atomB).name; nameC = GetAtomPara(atomC).name;
	biliA = (double)A / N;
	biliB = (double)B / N;
	biliC = (double)C / N;
	
	R = calloc(N*N,sizeof(double));
	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);


	Distance(x,y,z,R,N);
	for(i=0;i<N;i++){
		x[i] = x[i]*a0/(Rmin*sqrt(2));
		y[i] = y[i]*a0/(Rmin*sqrt(2));
		z[i] = z[i]*a0/(Rmin*sqrt(2));
	}
	Distance(x,y,z,R,N);

	SetEnergyPow(atomA,atomB,atomC);
	setupJohnson();
//	setEnergyNewModel(atomA,atomB,R,N);

	//初始化
	for(i=0;i<POPSIZE;i++){
		pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomC,N,a0);
	//	pop[i].energy = QSCEnergyNewModel(pop[i].chrom,atomA,atomB,N);
	}
	best_pop.chrom = calloc(N,sizeof(int));
	worst_pop.chrom = calloc(N,sizeof(int));
	changeIndividual(&best_pop,&(pop[0]),N);
	changeIndividual(&worst_pop,&(pop[0]),N);
	for(i=1;i<POPSIZE;i++)
	{
		if(pop[i].energy<best_pop.energy)
			changeIndividual(&best_pop,&(pop[i]),N);
		if(pop[i].energy>worst_pop.energy)
			changeIndividual(&worst_pop,&(pop[i]),N);
	}

	E0 = best_pop.energy;
	E1 = E0;
	
	printf("\nInit GA...\n");
	printf("shape=%s\tN=%d\n",shape,N);
	printf("POPSIZE=%d\tpc=%lf\tpm=%lf\n",POPSIZE,pc,pm);
	printf("biliA=%.3f\tbiliB=%.3f\tbiliC=%.3f\n",biliA,biliB,biliC);
	printf("%s=%d\t%s=%d\t%s=%d\n",nameA,A,nameB,B,nameC,C);	
	printf("best=%lf\tworst=%lf\n",best_pop.energy,worst_pop.energy);


	printf("\nStart GA...\n");
	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\energy.txt");
	fp = fopen(Line_Date,"w");
	fprintf(fp,"\t%d\t%3.1lf\n",N,biliA);
	fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);
	while(1){
		E0 = best_pop.energy;
		clocks++;
	
		//选择
		select_operator(pop,N,POPSIZE,best_pop.energy,worst_pop.energy);
		//交叉
		crossover_operator(pop,N,POPSIZE,pc);
		//变异
		mutation_operator(pop,N,POPSIZE,pm);
		

		for(i=0;i<POPSIZE;i++)
			pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomC,N,a0);
		//	pop[i].energy = QSCEnergyNewModel(pop[i].chrom,atomA,atomB,N);
		
		adjustment_operator(rate,pop,R,atomA,atomB,atomC,energy,N,POPSIZE);
		
		//保留最优个体
		changeIndividual(&pop[0],&best_pop,N);

		changeIndividual(&best_pop,&(pop[0]),N);
		changeIndividual(&worst_pop,&(pop[0]),N);
		for(i=0;i<POPSIZE;i++)
		{
			if(pop[i].energy<best_pop.energy)
				changeIndividual(&best_pop,&(pop[i]),N);
			if(pop[i].energy>worst_pop.energy)
				changeIndividual(&worst_pop,&(pop[i]),N);
		}
		
		E1 = best_pop.energy;
		if( E1-E0<0.001 && E1-E0>-0.001){
			tig = 0;
			step ++;
		}	else {
			tig = 1;
			step = 0;
			IJKL++;
		}
		if(tig == 1 || (step!=0 && step%100 == 0))
		{	
			E0 = E1;
			printf("clocks=%d\tIJKL=%d\tstep=%d\n",clocks,IJKL,step);
			printf("  E   =%f\n",E0);
			fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);

			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\result.txt");			
			printResult3(best_pop.chrom,x,y,z,N,0,Line_Date);

			strcpy(Line_Date,Line_End);
			strcat(Line_Date,"\\Diamond.txt");
			printDiamond3(best_pop.chrom,x,y,z,atomA,atomB,atomC,N,Line_Date);
		}

		if(step>1000)
			break;

	}
	fclose(fp);
	free(R);
	free(worst_pop.chrom);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");			
	printResult3(best_pop.chrom,x,y,z,N,0,Line_Date);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\Diamond.txt");
	printDiamond3(best_pop.chrom,x,y,z,atomA,atomB,atomC,N,Line_Date);

	strcpy(Line_Date,Line_End);
	strcat(Line_Date,"\\result.txt");
	printData(Line_Date,Line_End,N);

	free(best_pop.chrom);
//	freeSetEnergyNewModel();
}

//选择算子
void select_operator(struct individual *pop,int N,int POPSIZE,double bestEnergy,double worstEnergy){
	int i,index;
	double p,sum=0;
	double *cvalue;
	struct individual *new_pop;
	double rou = 3;
	double a1=0;
	double a2 = 1;
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

//交叉算子
void crossover_operator(struct individual *pop,int N,int POPSIZE,double pc){
	int i,j;
	int length;
	int *index;
	int dNum[3];
	int point,point1,point2,point_left,point_right,temp;
	double p;
	int ch;
	
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

//变异算子
void mutation_operator(struct individual *pop,int N,int POPSIZE,double pm){
	int i,j;
	int point;
	int ch;
	double p;

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

void adjustment_operator(double rate,  struct individual *pop, double *R, ATOM atomA, ATOM atomB, ATOM atomC, PE energy, int N, int POPSIZE)
{
	int i,randN1,randN2,tempChrom;
	double tempE,r;
	struct individual *pbest;

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
			tempE = GetCutEnergyFunction(energy)(pbest->chrom,R,atomA,atomB,atomB,N,getLatticeParameter3(atomA,atomB,atomC)); 
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


//交换个体间的数据
void changeIndividual(struct individual *one,struct individual *two,int N){
	int i;
	one->energy=two->energy;
	one->number=two->number;
	one->value=two->value;
	for(i=0;i<N;i++){
		one->chrom[i] = two->chrom[i];
	}
}

void GA_InitWithMixing(char *shape,int N, ATOMNUM atomNum,ALLOY alloy, PE energy, GAPARA para,char *Output)
{
	int i;
	COOD cood;
	struct individual *pop;

	cood.x = calloc(N,sizeof(double));
	cood.y = calloc(N,sizeof(double));
	cood.z = calloc(N,sizeof(double));
	cood.N = N;

	ReadCood(shape,N,cood.x,cood.y,cood.z);
	orderCoodFromCore(cood.x,cood.y,cood.z,N);

	pop = calloc(para.popSize,sizeof(INDIVIDUAL));

	srand((unsigned)time(NULL));
	for(i=0;i<para.popSize;i++){
		pop[i].number = i;
		pop[i].chrom = calloc(N,sizeof(int));
		MixNoteInt(pop[i].chrom, N, atomNum);
	}

	GA_Start(shape,N,atomNum,pop,cood,alloy,energy,para,Output);

	free(cood.x);
	free(cood.y);
	free(cood.z);
	for(i = 0;i < para.popSize;i++)
		free(pop[i].chrom);
	free(pop);
	
}

void GA_Start(char *shape, int N, ATOMNUM atomNum, INDIVIDUAL *pop, COOD cood, ALLOY alloy, PE energy, GAPARA para, char *output)
{
	int i;
	FILE *fp;
	double E0,E1;
	int clocks = 0, IJKL = 0, tig = 0, step = 0;
	char *Line_End;
	char Line_Date[100];
	struct individual best_pop;
	struct individual worst_pop;
	double a0;
	double biliA,biliB,biliC;
	char *nameA,*nameB,*nameC;
	COODDIS dis;

	srand((unsigned)time(NULL));
	dis.R = calloc(N*N,sizeof(double));
	a0 = getLatticeParameter(alloy.atoms,alloy.atomTypeCount);

	//	Line_End = StoragePath(shape,N,A,B,atomA,atomB,atomC,Output);


//	Distance1(cood,&dis);
//	for(i=0;i<N;i++){
//		x[i] = x[i]*a0/(Rmin*sqrt(2));
//		y[i] = y[i]*a0/(Rmin*sqrt(2));
//		z[i] = z[i]*a0/(Rmin*sqrt(2));
//	}
//	Distance(x,y,z,R,N);
//
//	SetEnergyPow(atomA,atomB,atomC);
//	setupJohnson();

//
//	//初始化
//	for(i=0;i<POPSIZE;i++){
//		pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomC,N,a0);
//	//	pop[i].energy = QSCEnergyNewModel(pop[i].chrom,atomA,atomB,N);
//	}
//	best_pop.chrom = calloc(N,sizeof(int));
//	worst_pop.chrom = calloc(N,sizeof(int));
//	changeIndividual(&best_pop,&(pop[0]),N);
//	changeIndividual(&worst_pop,&(pop[0]),N);
//	for(i=1;i<POPSIZE;i++)
//	{
//		if(pop[i].energy<best_pop.energy)
//			changeIndividual(&best_pop,&(pop[i]),N);
//		if(pop[i].energy>worst_pop.energy)
//			changeIndividual(&worst_pop,&(pop[i]),N);
//	}
//
//	E0 = best_pop.energy;
//	E1 = E0;
//	
//	printf("\nInit GA...\n");
//	printf("shape=%s\tN=%d\n",shape,N);
//	printf("POPSIZE=%d\tpc=%lf\tpm=%lf\n",POPSIZE,pc,pm);
//	printf("biliA=%.3f\tbiliB=%.3f\tbiliC=%.3f\n",biliA,biliB,biliC);
//	printf("%s=%d\t%s=%d\t%s=%d\n",nameA,A,nameB,B,nameC,C);	
//	printf("best=%lf\tworst=%lf\n",best_pop.energy,worst_pop.energy);
//
//
//	printf("\nStart GA...\n");
//	strcpy(Line_Date,Line_End);
//	strcat(Line_Date,"\\energy.txt");
//	fp = fopen(Line_Date,"w");
//	fprintf(fp,"\t%d\t%3.1lf\n",N,biliA);
//	fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);
//	while(1){
//		E0 = best_pop.energy;
//		clocks++;
//	
//		//选择
//		select_operator(pop,N,POPSIZE,best_pop.energy,worst_pop.energy);
//		//交叉
//		crossover_operator(pop,N,POPSIZE,pc);
//		//变异
//		mutation_operator(pop,N,POPSIZE,pm);
//		
//
//		for(i=0;i<POPSIZE;i++)
//			pop[i].energy = GetCutEnergyFunction(energy)(pop[i].chrom,R,atomA,atomB,atomC,N,a0);
//		//	pop[i].energy = QSCEnergyNewModel(pop[i].chrom,atomA,atomB,N);
//		
//		adjustment_operator(rate,pop,R,atomA,atomB,atomC,energy,N,POPSIZE);
//		
//		//保留最优个体
//		changeIndividual(&pop[0],&best_pop,N);
//
//		changeIndividual(&best_pop,&(pop[0]),N);
//		changeIndividual(&worst_pop,&(pop[0]),N);
//		for(i=0;i<POPSIZE;i++)
//		{
//			if(pop[i].energy<best_pop.energy)
//				changeIndividual(&best_pop,&(pop[i]),N);
//			if(pop[i].energy>worst_pop.energy)
//				changeIndividual(&worst_pop,&(pop[i]),N);
//		}
//		
//		E1 = best_pop.energy;
//		if( E1-E0<0.001 && E1-E0>-0.001){
//			tig = 0;
//			step ++;
//		}	else {
//			tig = 1;
//			step = 0;
//			IJKL++;
//		}
//		if(tig == 1 || (step!=0 && step%100 == 0))
//		{	
//			E0 = E1;
//			printf("clocks=%d\tIJKL=%d\tstep=%d\n",clocks,IJKL,step);
//			printf("  E   =%f\n",E0);
//			fprintf(fp,"%6d\t%6d\t%lf\n",clocks,IJKL,E0);
//
//			strcpy(Line_Date,Line_End);
//			strcat(Line_Date,"\\result.txt");			
//			printResult3(best_pop.chrom,x,y,z,N,0,Line_Date);
//
//			strcpy(Line_Date,Line_End);
//			strcat(Line_Date,"\\Diamond.txt");
//			printDiamond3(best_pop.chrom,x,y,z,atomA,atomB,atomC,N,Line_Date);
//		}
//
//		if(step>1000)
//			break;
//
//	}
//	fclose(fp);
//	free(R);
//	free(worst_pop.chrom);
//
//	strcpy(Line_Date,Line_End);
//	strcat(Line_Date,"\\result.txt");			
//	printResult3(best_pop.chrom,x,y,z,N,0,Line_Date);
//
//	strcpy(Line_Date,Line_End);
//	strcat(Line_Date,"\\Diamond.txt");
//	printDiamond3(best_pop.chrom,x,y,z,atomA,atomB,atomC,N,Line_Date);
//
//	strcpy(Line_Date,Line_End);
//	strcat(Line_Date,"\\result.txt");
//	printData(Line_Date,Line_End,N);
//
//	free(best_pop.chrom);
}