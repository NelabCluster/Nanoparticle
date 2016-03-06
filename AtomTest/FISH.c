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
	
	Distance1(&instance->cood,&instance->dis);
	for(i=0;i<instance->N;i++){
		instance->cood.x[i] = instance->cood.x[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.y[i] = instance->cood.y[i]*a0/(instance->dis.Rmin*sqrt(2));
		instance->cood.z[i] = instance->cood.z[i]*a0/(instance->dis.Rmin*sqrt(2));
	}
	Distance1(&instance->cood,&instance->dis);
	
	// 计算种群中每个个体的能量
	for( i = 0; i < instance->para.popSize; i++)
		instance->pop[i].energy = GetEnergyFunction1(instance->energyType)(instance->pop[i].chrom,instance->dis.R,instance->alloy.atoms,instance->alloy.atomTypeCount,instance->N);
	
	FISH_RefreshBest( instance );
	
	if( instance->para.selfAdaption == 1 )
		FISH_SelAdaption( instance );

	FISH_RefreshVisual( instance );
	
	// 打印鱼群算法实例的基本信息
	printf( "\nInit FISH...\n" );
	FISH_PrintMsg( instance );

	printf( "\nStart FISH...\n" );
	FISH_EnergyFile( instance, &fp, Line_End );
	while(1)
	{
		instance->clocks ++;
		
		if( instance->para.selfAdaption == 1 )
			FISH_SelAdaption( instance );

		for( i = 0; i < instance->para.popSize; i ++ )
		{
			state = FISH_Follow( instance, i );
			if( state != 1)
			{
				
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

void FISH_Prey( FISHInstance *instance, int i )
{ }

void FISH_Swarm( FISHInstance *instance, int i )
{ }

void FISH_SelAdaption( FISHInstance *instance )
{ }

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
{ }

void FISH_EnergyFile(FISHInstance *instance, FILE **fp, char *output)
{ }

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