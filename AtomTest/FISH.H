#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"

/***************************************
* @Name: 鱼群算法参数
****************************************/
typedef struct
{
	int popSize;				//种群个数
	double delta;				//拥挤度因子
	char selfAdaption;			//是否需要参数自适应；1：需要，0：不需要
	int tryNumber;				//觅食尝试次数，当selfAdaption是1时无效
	int visual;					//可视距离，当selfAdaption是1时无效
	int generation;				//总迭代次数
	char needOrderCood;			//是否需要排序；1：需要，0：不需要
	char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
	double adjustmentRate;		//调整概率，当needAdjustment是1时有效
} FISHPARA;

/***************************************
* @Name: 鱼群算法默认参数值
****************************************/
static const FISHPARA _defFISHPara = {
	100,0.8,1,100,100,2000,0,0,1.0
};


/***************************************
* @Name: 鱼群算法个体
****************************************/
struct FISHIndividual
{	
	int *chrom;							//原子排布序列		
	double energy;						//能量
};
typedef struct FISHIndividual FISHINDIVIDUAL;

/***************************************
* @Name: 鱼群算法种群可视表
****************************************/
typedef struct
{
	char *table;
	int visual;
} FISVisualTable;

/***************************************
* @Name: 鱼群算法实例
****************************************/
typedef struct
{
	char *shape;				//构型名称
	int N;						//尺寸
	PE	energyType;				//势能类型
	ALLOY alloy;				//合金类型
	ATOMNUM atomNum;			//各类金属的原子个数
	COOD cood;					//坐标值
	COODDIS dis;				//原子间距离
	FISHPARA para;				//鱼群算法参数
	FISVisualTable visualTable; //可视表
	FISHINDIVIDUAL *pop;		//种群
	FISHINDIVIDUAL best;		//最优个体
	int clocks;					//迭代次数
	int IJKL;					//找到最优值的次数
} FISHInstance;

/***************************************
* @Name:FISHPara_Init
* @Purpose：给para赋值默认的值
			默认值为静态变量 _defFISHPara
* @param：FISHPARA *para --- 鱼群算法参数
* @return：void --- 
****************************************/
void FISHPara_Init(FISHPARA *para);

/***************************************
* @Name:FISH2_InitWithMixing
* @Purpose：鱼群算法 二合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            PE energy --- 势能类型
            FISHPARA *para --- 鱼群算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void FISH2_InitWithMixing(char *shape,int N,int A,
						ATOM atomA,ATOM atomB,
						PE energy,FISHPARA *para,char *output);

/***************************************
* @Name:FISH3_InitWithMixing
* @Purpose：鱼群算法 三合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能类型
			FISHPARA *para --- 鱼群算法参数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void FISH3_InitWithMixing(char *shape,int N,int A,int B,
						ATOM atomA,ATOM atomB,ATOM atomC,
						PE energy,FISHPARA *para,char *output);

/***************************************
* @Name:FISH_InitWithMixing
* @Purpose：鱼群算法 多合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 原子总个数
            ATOMNUM *atomNum --- 各原子数量 
            ALLOY *alloy --- 合金类型
            PE energy --- 能量类型
            FISHPARA *para --- 鱼群算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void FISH_InitWithMixing(char *shape, int N, 
					   ATOMNUM *atomNum, ALLOY *alloy, 
					   PE energy, FISHPARA *para,char *output);

/***************************************
* @Name:FISH_Start
* @Purpose：运行鱼群算法实例，会在[output]文件夹下生成结果文件和分析文件
* @param：FISHInstance *instance --- 鱼群算法实例
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void FISH_Start( FISHInstance *instance, char *output);

char FISH_Follow( FISHInstance *instance, int index );

void FISH_Prey( FISHInstance *instance, int index );

void FISH_Swarm( FISHInstance *instance, int index );

void FISH_SelAdaption( FISHInstance *instance );

void FISH_RefreshBest( FISHInstance *instance );

void FISH_RefreshVisual( FISHInstance *instance );

void FISH_RefreshVisualAtIndex( FISHInstance *instance, int index );

/***************************************
* @Name:FISH_PrintMsg
* @Purpose：打印鱼群算法实例的基本信息，如构型名称，原子总数，原子类型，各类型个数和最优值和最差值
* @param：FISHInstance *instance --- 鱼群算法实例
* @return：void --- 
****************************************/
void FISH_PrintMsg(FISHInstance *instance);

/***************************************
* @Name:FISH_EnergyFile
* @Purpose：输出能量文件，记录每一代有效的最优值
			文件路径为：[output]\\energy.txt
			当输出结束之后，需要手动fclose(fp)
* @param：FISHInstance *instance --- 鱼群算法实例
            FILE **fp --- 文件指针的地址
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void FISH_EnergyFile(FISHInstance *instance, FILE **fp, char *output);


/***************************************
* @Name:FISH_ResultFile
* @Purpose：在在文件夹路径下生成结果文件:
			结果文件 - [output]\\result.txt
			绘图文件 - [output]\\Diamond.txt
* @param：FISHInstance *instance --- 鱼群算法实例
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void FISH_ResultFile( FISHInstance *instance, char *output );

/***************************************
* @Name:FISHIndividual_Init
* @Purpose：初始化个体实例，为其开辟所需的空间
			每一个个体实例都需要初始化
* @param：FISHINDIVIDUAL *one --- 鱼群算法个体实例
            int N --- 总原子数
* @return：void --- 
****************************************/
void FISHIndividual_Init( FISHINDIVIDUAL *one, int N );


/***************************************
* @Name:FISHIndividual_Copy
* @Purpose：将源个体实例from中的值赋值给目标个体实例
			to 和 from 都必须被初始化过，即带入过 FISHIndividual_Init
* @param：FISHINDIVIDUAL *to --- 目标个体实例
            FISHINDIVIDUAL *from --- 源个体实例
            int N --- 原子总数
* @return：void --- 
****************************************/
void FISHIndividual_Copy( FISHINDIVIDUAL *to, FISHINDIVIDUAL *from, int N );


/***************************************
* @Name:FISHIndividual_Free
* @Purpose：释放个体实例初始化中开辟的空间
			当一个个体实例one运行 FISHIndividual_Init(one) 后，当不使用时候必须执行 FISHIndividual_Free(one)
* @param：INDIVIDUAL *one --- 鱼群算法个体实例
* @return：void --- 
****************************************/
void FISHIndividual_Free( FISHINDIVIDUAL *one );

int FISHIndividual_Dist( FISHINDIVIDUAL *one, FISHINDIVIDUAL *two, int N );

void FISVisualTable_Init( FISVisualTable *table, int popSize );

void FISVisualTable_Free( FISVisualTable *table );