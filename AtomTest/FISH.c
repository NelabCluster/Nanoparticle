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
	
	//��ʼ����Ⱥ�㷨ʵ��
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

	//���ݹ��ͺ�ԭ������ȡ����
	ReadCood1(shape,N,&instance.cood);
	
	//�����Ҫ������������ﵽ����������
	if(instance.para.needOrderCood == 1)
		orderCoodFromCore( instance.cood.x, instance.cood.y, instance.cood.z, instance.N );
	
	//����������ӣ���Ȼÿ���������ֵ��һ����ÿ���������������ʼ��
	srand((unsigned)time(NULL));
	for( i = 0; i < instance.para.popSize;i++)
		MixNoteInt( instance.pop[i].chrom, instance.N, &instance.atomNum );

	//��ʼ������Ⱥ�㷨
	FISH_Start( &instance, output );
	
	//�ͷ���Ⱥ�㷨ʵ���п��ٵĿռ�
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
	
	Distance1(&instance->cood,&instance->dis);
	for(i=0;i<instance->N;i++){
		instance->cood.x[i] = instance->cood.x[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.y[i] = instance->cood.y[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.z[i] = instance->cood.z[i]*a0/(instance->dis.Rmin*sqrt(2));
	}
	Distance1(&instance->cood,&instance->dis);
	
	// ������Ⱥ��ÿ�����������
	for( i = 0; i < instance->para.popSize; i++)
		instance->pop[i].energy = GetCutEnergyFunction1(instance->energyType)(instance->pop[i].chrom,instance->dis.R,&instance->alloy,instance->N);
	
	FISH_RefreshBest( instance );
		
	E0 = instance->best.energy;
	E1 = E0;

	if( instance->para.selfAdaption == 1 )
		FISH_SelAdaption( instance );

	FISH_RefreshVisual( instance );
	
	// ��ӡ��Ⱥ�㷨ʵ���Ļ�����Ϣ
	printf( "\nInit FISH...\n" );
	FISH_PrintMsg( instance );

	printf( "\nStart FISH...\n" );
	FISH_EnergyFile( instance, &fp, Line_End );
	while(1)
	{
		instance->clocks ++;
		
		if( instance->para.selfAdaption == 1 )
			FISH_SelAdaption( instance );

		FISH_RefreshVisual( instance );
		
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

		if( instance->clocks >= instance->para.generation ) break; 
	}

	free( instance->dis.R );
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

		if( visualSize / popSize < instance->para.delta )
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
		
		tempE = GetCutEnergyFunction1(instance->energyType)(one->chrom,instance->dis.R,&instance->alloy,N); 
		
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

char FISH_Swarm( FISHInstance *instance, int i )
{ return 0; }

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

void FISH_PrintMsg(FISHInstance *instance)
{ 
	int i,tempNum;
	ATOM tempATOM;

	printf( "shape=%s\tN=%d\n", instance->shape, instance->N );
	printf( "POPSIZE=%d\ttryNumber=%d\tvisual=%lf\n", instance->para.popSize, instance->para.tryNumber, instance->para.visual );
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf( "%s=%.3f\t", GetAtomPara(tempATOM).name, (double)tempNum / instance->N );
	}
	printf("\n");
	for( i = 0; i < instance->alloy.atomTypeCount; i++ )
	{
		tempATOM = instance->alloy.atoms[i];
		tempNum = instance->atomNum.numberOfAtom[i];
		printf( "%s=%d\t", GetAtomPara(tempATOM).name, tempNum );
	}
	printf( "\n" );
	printf( "best=%lf\n", instance->best.energy );
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
{ }

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