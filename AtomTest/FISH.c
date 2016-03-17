#include "FISH.h"

void FISHPara_Init(FISHPARA *para)
{
	memcpy( para, &_defFISHPara, sizeof(*para) );
}

void FISH2_InitWithMixing(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,FISHPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,END);
	AtomNum_Init(&atomNum,A,N-A,END);
	FISH_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);
}

void FISH3_InitWithMixing(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,FISHPARA *para,char *output)
{
	ALLOY alloy;
	ATOMNUM atomNum;

	Alloy_Init(&alloy,atomA,atomB,atomC,END);
	AtomNum_Init(&atomNum,A,B,N-A-B,END);
	FISH_InitWithMixing(shape,N,&atomNum,&alloy,energy,para,output);
	
	Alloy_Free(&alloy);
	AtomNum_Free(&atomNum);	
}

void FISH_InitWithMixing(char *shape, int N, ATOMNUM *atomNum, ALLOY *alloy, PE energy, FISHPARA *para,char *output)
{
	int i;
	FISHInstance instance;
	
	//初始化鱼群算法实例
	instance.shape = shape;
	instance.N = N;
	instance.energyType = energy;
	Alloy_Copy( &instance.alloy, alloy );
	AtomNum_Copy( &instance.atomNum, atomNum );
	Cood_Init( &instance.cood, N );
	instance.para = *para;
	FISVisualTable_Init( &instance.visualTable, para->popSize );
	instance.pop = calloc( para->popSize, sizeof(FISHINDIVIDUAL) );
	for( i = 0; i < para->popSize; i++)
		FISHIndividual_Init( &instance.pop[i], N );
	FISHIndividual_Init( &instance.best, N );
	instance.clocks = 0;
	instance.IJKL = 0;

	//根据构型和原子数读取坐标
	ReadCood1(shape,N,&instance.cood);
	
	//如果需要排序坐标则从里到外排序坐标
	if(instance.para.needOrderCood == 1)
		orderCoodFromCore( instance.cood.x, instance.cood.y, instance.cood.z, instance.N );
	
	//设置随机种子，不然每次启动随机值都一样，每个个体随机产生初始解
	srand((unsigned)time(NULL));
	for( i = 0; i < instance.para.popSize;i++)
		MixNoteInt( instance.pop[i].chrom, instance.N, &instance.atomNum );

	//开始进行鱼群算法
	FISH_Start( &instance, output );
	
	//释放鱼群算法实例中开辟的空间
	FISHIndividual_Free( &instance.best );
	for(i = 0;i < instance.para.popSize;i++)
		FISHIndividual_Free( &instance.pop[i] );
	free( instance.pop );
	FISVisualTable_Free( &instance.visualTable );
	Cood_Free( &instance.cood );
	AtomNum_Free( &instance.atomNum );
	Alloy_Free( &instance.alloy );
}

void FISH_Start(FISHInstance *instance, char *output)
{ 
	FILE *fp = NULL;
	double a0,E0,E1;
	int i, tig = 0, step = 0;
	char state;
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
	
	FISH_RefreshBest( instance );
		
	E0 = instance->best.energy;
	E1 = E0;

	if( instance->para.selfAdaption == 1 )
		FISH_SelAdaption( instance );

	FISH_RefreshVisual( instance );
	
	// 打印鱼群算法实例的基本信息
	printf( "\nInit FISH...\n" );
	FISH_PrintMsg( instance, Line_End );

	printf( "\nStart FISH...\n" );
	FISH_EnergyFile( instance, &fp, Line_End );
	while(1)
	{
		instance->clocks ++;
		
		if( instance->para.selfAdaption == 1 )
			FISH_SelAdaption( instance );

		FISH_RefreshVisual( instance );
		//	printf("try=%d\n",instance->para.tryNumber);
		//	printf("1\n");
		for( i = 0; i < instance->para.popSize; i ++ )
		{
			state = FISH_Follow( instance, i );
			if( state != 1)
			{
				state = FISH_Prey( instance, i );
				
				if( state != 1 )
					state = FISH_Swarm( instance, i );
			}
		}
		//	printf("2,step=%d\n",step);
		
		FISH_RefreshBest( instance );
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

			FISH_EnergyFile(instance,&fp,Line_End);
			FISH_ResultFile(instance,Line_End);	
		}



		if( step>instance->para.convergenceGenerations ||  instance->clocks >= instance->para.generation )
			break; 
	}



	free( instance->dis.R );
	Energy_Free( instance->energyType );

	printf("\nEnd GA...\n");
	FISH_PrintMsg( instance, Line_End );
	FISH_ResultFile( instance, Line_End );

	printData( Line_End, instance->N );
	free( Line_End );
}

char FISH_Follow( FISHInstance *instance, int index )
{
	int i,popSize,visualSize,bestInVisualIndex;
	
	popSize = instance->para.popSize;
	bestInVisualIndex = -1;
	for( i = 0; i < popSize; i++ )
	{
		if( i == index ) continue;

		if( instance->visualTable.table[ index*popSize+i ] == 1 )
		{
			if( bestInVisualIndex == -1 || instance->pop[i].energy <= instance->pop[bestInVisualIndex].energy )
				bestInVisualIndex = i;
		}
	}
	
	if( bestInVisualIndex != -1 && instance->pop[bestInVisualIndex].energy >= instance->pop[index].energy )
		bestInVisualIndex = -1;
	
	if( bestInVisualIndex != -1 )
	{
		visualSize = 0;
		for( i = 0; i < popSize; i++ )
		{
			if( i == bestInVisualIndex ) continue;

			if( instance->visualTable.table[ bestInVisualIndex*popSize+i ] == 1 )
				visualSize ++;
		}

		if( visualSize*1.0 / popSize < instance->para.delta )
		{
			FISHIndividual_Copy( instance->pop+index, instance->pop+bestInVisualIndex, instance->N );
			FISH_RefreshVisualAtIndex( instance, index );
		}
		else
			bestInVisualIndex = -1;
	}

	return ( bestInVisualIndex == -1 )?0:1;
}

char FISH_Prey( FISHInstance *instance, int index )
{ 
	int i,randN1,randN2,tempChrom,POPSIZE,N;
	double bestEnergyAfterSwap,tempE;
	int swap[2];
	FISHINDIVIDUAL *one;

	POPSIZE = instance->para.popSize;
	N = instance->N;
	one = instance->pop + index;
	bestEnergyAfterSwap = 0;
	swap[0] = 0; swap[1] = 0;
	
	for( i = 0; i < instance->para.tryNumber; i++ )
	{
		randN1 = RANDINT( N );
		randN2 = RANDINT( N );
		while( one->chrom[randN1] == one->chrom[randN2] )
		{
			randN1 = RANDINT( N );
			randN2 = RANDINT( N );
		}
		tempChrom = one->chrom[randN1];
		one->chrom[randN1] = one->chrom[randN2];
		one->chrom[randN2] = tempChrom;
		
		tempE = GetCutEnergyFunction(instance->energyType)(one->chrom,instance->dis.R,&instance->alloy,N); 
		
		if( tempE < bestEnergyAfterSwap )
		{
			bestEnergyAfterSwap = tempE;
			swap[0] = randN1;
			swap[1] = randN2;
		}

		tempChrom = one->chrom[randN1];
		one->chrom[randN1] = one->chrom[randN2];
		one->chrom[randN2] = tempChrom;
	}

	if( bestEnergyAfterSwap  < one->energy )
	{
		one->energy = bestEnergyAfterSwap;
		tempChrom = one->chrom[swap[0]];
		one->chrom[swap[0]] = one->chrom[swap[1]];
		one->chrom[swap[1]] = tempChrom;
		FISH_RefreshVisualAtIndex( instance, index );
	} else
	{
		bestEnergyAfterSwap = 0;
		swap[0] = 0; swap[1] = 0;
	}
	
	return ( bestEnergyAfterSwap == 0 )?0:1;
}

char FISH_Swarm( FISHInstance *instance, int index )
{
	int i,j,N,popSize,nearNum,atomTypeCount,maxIndex;
	FISHINDIVIDUAL center;
	int *numOfAtoms;
	
	FISHIndividual_Init( &center, instance->N );
	N = instance->N;
	popSize = instance->para.popSize;
	numOfAtoms = calloc( instance->alloy.atomTypeCount, sizeof(int) );
	atomTypeCount = instance->alloy.atomTypeCount;

	nearNum = 0;
	for( i = 0; i < popSize; i++ )
	{
		if( i == index ) continue;
		if( instance->visualTable.table[index*popSize + i] == 1 )
			nearNum ++;
	}

	if( nearNum == 0 )
		return 0;

	if( nearNum*1.0 / popSize >= instance->para.delta )
		return 0;
	
	for( i = 0; i < N; i++ )
	{
		memset( numOfAtoms, 0, atomTypeCount*sizeof(int) );
		for( j = 0; j < popSize; j++ )
		{
			if( j == index ) continue;
			if( instance->visualTable.table[index*popSize + j] == 1 )
			{
				numOfAtoms[ instance->pop[j].chrom[i] ] ++;
			}
		}
		
		maxIndex = 0;
		for( j = 1; j < atomTypeCount; j++ )
		{
			if( numOfAtoms[j] > numOfAtoms[maxIndex] )
				maxIndex = j;
		}
		center.chrom[i] = maxIndex;
	}
	
	//	FISHIndividual_Print( &center, &instance->atomNum, N );
	FISHIndividual_Adjust( &center, &instance->atomNum, N );
	//	FISHIndividual_Print( &center, &instance->atomNum, N );

	center.energy = GetCutEnergyFunction(instance->energyType)(center.chrom,instance->dis.R,&instance->alloy,N); 
	
	if( center.energy < instance->pop[index].energy )
	{
		FISHIndividual_Copy( instance->pop+index, &center, instance->N );
		FISH_RefreshVisualAtIndex( instance, index );
	}

	free( numOfAtoms );
	FISHIndividual_Free( &center );
	return 0; 
}

void FISH_SelAdaption( FISHInstance *instance )
{
	int count,bb,tt,k;

	if( instance->clocks < 300 )
		instance->para.tryNumber = 10;
	else
		instance->para.tryNumber = (int)( 100*instance->clocks / instance->para.generation );
			
	if( instance->clocks == 0 )
		instance->para.visual = (int)( instance->N / 8 );

	if( instance->clocks%200 == 0 )
	{
		count = 0;
		bb = RANDINT( instance->para.popSize );
		for( k = 0; k < instance->para.popSize; k++ ){
			tt = FISHIndividual_Dist( instance->pop+bb, instance->pop+k, instance->N );
			count = count + tt;
		}
		instance->para.visual = (int) ( 0.5*count / instance->para.popSize );
	}
	//	instance->para.visual = 30;
}

void FISH_RefreshBest( FISHInstance *instance )
{ 
	int i;

	for( i = 0; i < instance->para.popSize; i++ )
		if( instance->pop[i].energy <= instance->best.energy )
			FISHIndividual_Copy( &instance->best, instance->pop+i, instance->N );
}

void FISH_RefreshVisual( FISHInstance *instance )
{
	int i,j,visual,popSize,dis;

	if( instance->visualTable.visual == instance->para.visual ) return;

	visual = instance->para.visual;
	popSize = instance->para.popSize;
	for( i =0; i < popSize; i++ )
		for( j = i+1; j < popSize; j++ )
		{
			dis = FISHIndividual_Dist( instance->pop+i, instance->pop+j, instance->N );
			if( dis < visual )
			{
				instance->visualTable.table[ i*popSize + j ] = 1;
				instance->visualTable.table[ j*popSize + i ] = 1;
			}
		}
}

void FISH_RefreshVisualAtIndex( FISHInstance *instance, int index )
{ 
	int i,popSize,dis;
	
	popSize = instance->para.popSize;
	for( i = 0; i < popSize; i++ )
	{
		if( i == index ) continue;

		dis = FISHIndividual_Dist( instance->pop+index, instance->pop+i, instance->N );
		if( dis < instance->visualTable.visual)
		{
			instance->visualTable.table[ index*popSize+i ] = 1;
			instance->visualTable.table[ i*popSize+index ] = 1;
		} else
		{
			instance->visualTable.table[ index*popSize+i ] = 0;
			instance->visualTable.table[ i*popSize+index ] = 0;			
		}
	}
}

void FISH_PrintMsg(FISHInstance *instance,  char *output )
{ 
	int i,tempNum;
	ATOM tempATOM;
	FILE *fp;
	char Line_Date[200];

	strcpy( Line_Date, output );
	strcat( Line_Date, "\\FISH.txt" );

	fp = fopen( Line_Date, "w" );

	printf( "shape=%s\tN=%d\n", instance->shape, instance->N );
	printf( "POPSIZE=%d\ttryNumber=%d\tvisual=%lf\n", instance->para.popSize, instance->para.tryNumber, instance->para.visual );
	fprintf( fp, "shape=%s\tN=%d\n", instance->shape, instance->N );
	fprintf( fp, "POPSIZE=%d\ttryNumber=%d\tvisual=%lf\n", instance->para.popSize, instance->para.tryNumber, instance->para.visual );

	
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf( "%s=%.3f\t", GetAtomPara(tempATOM).name, (double)tempNum / instance->N );
		fprintf( fp, "%s=%.3f\t", GetAtomPara(tempATOM).name, (double)tempNum / instance->N );
	}
	printf( "\n" );
	fprintf( fp, "\n" );
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf( "%s=%d\t", GetAtomPara(tempATOM).name, tempNum );
		fprintf( fp, "%s=%d\t", GetAtomPara(tempATOM).name, tempNum );
	}
	printf( "\n" );
	printf( "best=%lf\n", instance->best.energy );
	fprintf( fp, "\n" );
	fprintf( fp, "best=%lf\n", instance->best.energy );
}

void FISH_EnergyFile(FISHInstance *instance, FILE **fp, char *output)
{
 	char Line_Date[200];

	if ( *fp == NULL )
	{
		strcpy( Line_Date, output );
		strcat( Line_Date, "\\energy.txt" );
		*fp = fopen( Line_Date, "w" );
		fprintf( *fp, "\t%d\n", instance->N );
	}
	fprintf( *fp, "%6d\t%6d\t%lf\n", instance->clocks, instance->IJKL, instance->best.energy );
}

void FISH_ResultFile( FISHInstance *instance, char *output )
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

void FISHIndividual_Init( FISHINDIVIDUAL *one, int N )
{ 
	one->chrom = calloc( N, sizeof(int) );
	one->energy = 0;
}

void FISHIndividual_Copy( FISHINDIVIDUAL *to, FISHINDIVIDUAL *from, int N )
{ 
	int i;

	to->energy = from->energy;
	for( i =0 ; i < N; i++ )
		to->chrom[i] = from->chrom[i];
}

void FISHIndividual_Free( FISHINDIVIDUAL *one )
{
	one->energy = 0;
	free( one->chrom );
	one->chrom = NULL;
}

int FISHIndividual_Dist( FISHINDIVIDUAL *one, FISHINDIVIDUAL *two, int N )
{
	int i, dis = 0;

	for( i = 0; i < N; i++ )
	{
		if( one->chrom[i] != two->chrom[i] ) 
			dis ++;
	}

	return dis;
} 

void FISHIndividual_Adjust( FISHINDIVIDUAL *one, ATOMNUM *atomNum, int N )
{
	int i,atomTypeCount,randIndex0,randIndex,tempChorm;
	ATOMNUM atomNumOfOne;
	
	atomTypeCount = atomNum->atomTypeCount;
	AtomNum_Copy( &atomNumOfOne, atomNum );

	for( i = 0; i < atomTypeCount; i++ )
		atomNumOfOne.numberOfAtom[i] = 0;

	for( i = 0; i < N; i++)
		atomNumOfOne.numberOfAtom[one->chrom[i]] ++;
	
	for( i = 0; i < atomTypeCount; i++ )
	{
		while( atomNumOfOne.numberOfAtom[i] < atomNum->numberOfAtom[i] )
		{
			if( atomNumOfOne.numberOfAtom[i] < atomNum->numberOfAtom[i] )
			{
				randIndex = RANDINT( N );
				tempChorm = one->chrom[ randIndex ];
				while( atomNumOfOne.numberOfAtom[tempChorm] <= atomNum->numberOfAtom[tempChorm] )
				{
					randIndex = RANDINT( N );
					tempChorm = one->chrom[ randIndex ];	
				}
				one->chrom[ randIndex ] = i;
				atomNumOfOne.numberOfAtom[i] ++;
				atomNumOfOne.numberOfAtom[tempChorm]--;
			}
//			if( atomNumOfOne.numberOfAtom[i] > atomNum->numberOfAtom[i] )
//			{
//				randIndex0 = RANDINT( N );
//				while( one->chrom[ randIndex ] != i )
//					randIndex0 = RANDINT( N );
//
//				randIndex = RANDINT( N );
//				tempChorm = one->chrom[ randIndex ];
//				while( atomNumOfOne.numberOfAtom[tempChorm] >= atomNum->numberOfAtom[tempChorm] )
//				{
//					randIndex = RANDINT( N );
//					tempChorm = one->chrom[ randIndex ];	
//				}
//				one->chrom[ randIndex0 ] = tempChorm;
//				atomNumOfOne.numberOfAtom[i] --;
//				atomNumOfOne.numberOfAtom[tempChorm]++;				
//			}
		}
	}

	AtomNum_Free( &atomNumOfOne );
}

void FISHIndividual_Print( FISHINDIVIDUAL *one, ATOMNUM *atomNum, int N )
{
	int i,atomTypeCount;
	ATOMNUM atomNumOfOne;
	
	atomTypeCount = atomNum->atomTypeCount;
	AtomNum_Copy( &atomNumOfOne, atomNum );

	for( i = 0; i < atomTypeCount; i++ )
		atomNumOfOne.numberOfAtom[i] = 0;

	for( i = 0; i < N; i++)
		atomNumOfOne.numberOfAtom[one->chrom[i]] ++;
	
	for( i = 0; i < atomTypeCount; i++ )
		printf("%d\t",atomNumOfOne.numberOfAtom[i]);
	printf("\n");

	AtomNum_Free( &atomNumOfOne );
}

void FISVisualTable_Init( FISVisualTable *table, int popSize )
{
	table->table = calloc( popSize * popSize, sizeof(char) );
	table->visual = 0;
}

void FISVisualTable_Free( FISVisualTable *table )
{
	free( table->table );
	table->table = NULL;
	table->visual = 0;
}