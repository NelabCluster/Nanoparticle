#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"


/***************************************
* @Name: 遗传算法参数
****************************************/
typedef struct
{
	int popSize;				//种群个数
	double pc;					//交叉概率
	double pm;					//变异概率
	int convergenceGenerations;	//收敛代数
	char needOrderCood;			//是否需要排序；1：需要，0：不需要
	char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
	double adjustmentRate;		//调整概率，当needAdjustment是1时有效
} GAPARA;

/***************************************
* @Name: 遗传算法默认参数值
****************************************/
static const GAPARA _defGAPara = {
	100,0.8,0.0015,1000,1,1,1.0
};

/***************************************
* @Name: 遗传算法个体
****************************************/
struct individual
{	
	int number;			//保留，未使用
	int *chrom;			//原子排布序列		
	double energy;		//能量
	double value;		//适应值
};
typedef struct individual INDIVIDUAL;

/***************************************
* @Name: 遗传算法实例
****************************************/
typedef struct
{
	char *shape;		//构型名称
	int N;				//尺寸
	PE	energyType;		//势能类型
	ALLOY alloy;		//合金类型
	ATOMNUM atomNum;	//各类金属的原子个数
	COOD cood;			//坐标值
	COODDIS dis;		//原子间距离
	GAPARA para;		//遗传算法参数
	INDIVIDUAL *pop;	//种群
	INDIVIDUAL best;	//最优个体
	INDIVIDUAL worst;	//最差个体
	int clocks;			//迭代次数
	int IJKL;			//找到最优值的次数
} GAInstance;


/***************************************
* @Name:GAPara_Init
* @Purpose：给para赋值默认的值
			默认值为静态变量 _defGAPara
* @param：GAPARA *para --- 遗传算法参数
* @return：void --- 
****************************************/
void GAPara_Init(GAPARA *para);

/***************************************
* @Name:GA2_InitWithMixing
* @Purpose：遗传算法 二合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            PE energy --- 势能类型
            GAPARA *para --- 遗传算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA2_InitWithMixing(char *shape,int N,int A,
						ATOM atomA,ATOM atomB,
						PE energy,GAPARA *para,char *output);

/***************************************
* @Name:GA3_InitWithMixing
* @Purpose：遗传算法 三合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能类型
			GAPARA *para --- 遗传算法参数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA3_InitWithMixing(char *shape,int N,int A,int B,
						ATOM atomA,ATOM atomB,ATOM atomC,
						PE energy,GAPARA *para,char *output);

/***************************************
* @Name:GA_InitWithMixing
* @Purpose：遗传算法 多合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 原子总个数
            ATOMNUM *atomNum --- 各原子数量 
            ALLOY *alloy --- 合金类型
            PE energy --- 能量类型
            GAPARA *para --- 遗传算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA_InitWithMixing(char *shape, int N, 
					   ATOMNUM *atomNum, ALLOY *alloy, 
					   PE energy, GAPARA *para,char *output);

/***************************************
* @Name:GA_Start
* @Purpose：运行遗传算法实例，会在[output]文件夹下生成结果文件和分析文件
* @param：GAInstance *instance --- 遗传算法实例
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void GA_Start(GAInstance *instance, char *output);

/***************************************
* @Name:select_operator
* @Purpose：选择算子
			采用轮盘赌的方式选择出 popSize 个新个体 替换原有 instance->pop
* @param：GAInstance *instance --- 遗传算法实例
* @return：void --- 
****************************************/
void select_operator(GAInstance *instance);

/***************************************
* @Name:crossover_operator
* @Purpose：交叉算子
* @param：GAInstance *instance ---  遗传算法实例
* @return：void --- 
****************************************/
void crossover_operator(GAInstance *instance);

/***************************************
* @Name:mutation_operator
* @Purpose：变异算子
			随机交换两个不同类型的原子
* @param：GAInstance *instance --- 遗传算法实例 
* @return：void --- 
****************************************/
void mutation_operator(GAInstance *instance);

/***************************************
* @Name:adjustment_operator
* @Purpose：调整算子
			随机交换两个不同类型的原子，并计算交换后的能量，如果个体能量降低，则替换原个体
* @param：GAInstance *instance --- 遗传算法实例
* @return：void --- 
****************************************/
void adjustment_operator(GAInstance *instance);

/***************************************
* @Name:GA_FindBestAndWorst
* @Purpose：找到遗传算法实例中的最优个体和最差个体赋值给 instance->best 和 instance->worst
* @param：GAInstance *instance --- 遗传算法的实例
* @return：void --- 
****************************************/
void GA_FindBestAndWorst(GAInstance *instance);

/***************************************
* @Name:GA_PrintMsg
* @Purpose：打印遗传算法实例的基本信息，如构型名称，原子总数，原子类型，各类型个数和最优值和最差值
* @param：GAInstance *instance --- 遗传算法实例
* @return：void --- 
****************************************/
void GA_PrintMsg(GAInstance *instance);

/***************************************
* @Name:GA_EnergyFile
* @Purpose：输出能量文件，记录每一代有效的最优值
			文件路径为：[output]\\energy.txt
			当输出结束之后，需要手动fclose(fp)
* @param：GAInstance *instance --- 遗传算法实例
            FILE **fp --- 文件指针的地址
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void GA_EnergyFile(GAInstance *instance, FILE **fp, char *output);


/***************************************
* @Name:GA_ResultFile
* @Purpose：在在文件夹路径下生成结果文件:
			结果文件 - [output]\\result.txt
			绘图文件 - [output]\\Diamond.txt
* @param：GAInstance *instance --- 遗传算法实例
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void GA_ResultFile(GAInstance *instance, char *output);


/***************************************
* @Name:GAIndividual_Init
* @Purpose：初始化个体实例，为其开辟所需的空间
			每一个个体实例都需要初始化
* @param：INDIVIDUAL *one --- 遗传算法个体实例
            int N --- 总原子数
* @return：void --- 
****************************************/
void GAIndividual_Init(INDIVIDUAL *one,int N);


/***************************************
* @Name:changeIndividual
* @Purpose：将源个体实例from中的值赋值给目标个体实例
			to 和 from 都必须被初始化过，即带入过 GAIndividual_Init
* @param：struct individual *to --- 目标个体实例
            struct individual *from --- 源个体实例
            int N --- 原子总数
* @return：void --- 
****************************************/
void changeIndividual(struct individual *to,struct individual *from,int N);


/***************************************
* @Name:GAIndividual_Free
* @Purpose：释放个体实例初始化中开辟的空间
			当一个个体实例one运行 GAIndividual_Init(one) 后，当不使用时候必须执行 GAIndividual_Free(one)
* @param：INDIVIDUAL *one --- 遗传算法个体实例
* @return：void --- 
****************************************/
void GAIndividual_Free(INDIVIDUAL *one);


//以下未调试，不可使用
//初始核壳，壳层的层数为surfaceLayer 0:最外层
void GA3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);
void GA3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int POPSIZE,double pc,double pm,double rate,char *Output);
