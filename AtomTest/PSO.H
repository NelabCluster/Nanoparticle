#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"

/***************************************
* @Name: 粒子群算法参数
****************************************/
typedef struct
{
	int popSize;				//种群大小
	double w;					//惯性系数
	int maxGenerations;			//最大迭代次数
	int convergenceGenerations;	//收敛代数
	char needOrderCood;			//是否需要排序；1：需要，0：不需要
	char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
	double rate;				//调整概率，当needAdjustment是1时有效
} PSOPARA;

/***************************************
* @Name: 粒子群算法默认参数值
****************************************/
static const PSOPARA _defPSOPara = {
	100,0.25,5000,1000,1,1,1.0
};

/***************************************
* @Name: 交换子
****************************************/
struct PSOSwap
{
	int element[2];			
	struct PSOSwap *next;	
};
typedef struct PSOSwap PSOSWAP;

/***************************************
* @Name:粒子群个体
****************************************/
struct PSOIndividual
{
	int *chrom;						//解
	double energy;					//能量
	PSOSWAP *swapList;				//交换序
	struct PSOIndividual *pbest;	//个体最优
};
typedef struct PSOIndividual PSOINDIVIDUAL;

/***************************************
* @Name: 粒子群算法实例
****************************************/
typedef struct
{
	char *shape;			//构型名称
	int N;					//尺寸
	PE	energyType;			//势能类型
	ALLOY alloy;			//合金类型
	ATOMNUM atomNum;		//各类金属的原子个数
	COOD cood;				//坐标值
	COODDIS dis;			//原子间距离
	PSOPARA para;			//粒子群算法参数
	PSOINDIVIDUAL *pop;		//种群
	PSOINDIVIDUAL *gbest;	//全局最优个体
	int clocks;				//迭代次数
	int IJKL;				//找到最优值的次数
} PSOInstance;

/***************************************
* @Name:PSOPara_Init
* @Purpose：给para赋值默认的值
			默认值为静态变量 _defPSOPara
* @param：PSOPARA *para --- 粒子群算法参数
* @return：void --- 
****************************************/
void PSOPara_Init(PSOPARA *para);

/***************************************
* @Name:PSO2_InitWithMixing
* @Purpose：粒子群算法二合金初始结构随机分布
* @param：char *shape --- 形状
            int N --- 原子总数
            int A --- 原子A个数
            ATOM atomA --- 原子A
            ATOM atomB --- 原子B
			PE energy --- 势能类型
			PSOPARA *para --- 粒子群算法参数
            char *output --- 文件存储路径
* @return：void --- 
****************************************/
void PSO2_InitWithMixing(char *shape,int N,int A,
						 ATOM atomA,ATOM atomB,
						 PE energy,PSOPARA *para,char *output); 

/***************************************
* @Name:PSO3_InitWithMixing
* @Purpose：粒子群算法 三合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能类型
			GAPARA *para --- 粒子群算法参数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void PSO3_InitWithMixing( char *shape, int N, int A, int B, 
						 ATOM atomA, ATOM atomB, ATOM atomC, 
						 PE energy, PSOPARA *para, char *output );

/***************************************
* @Name:PSO_InitWithMixing
* @Purpose：粒子群算法 多合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 原子总个数
            ATOMNUM *atomNum --- 各原子数量 
            ALLOY *alloy --- 合金类型
            PE energy --- 能量类型
            PSOPARA *para --- 粒子群算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void PSO_InitWithMixing(char *shape,int N,
						ATOMNUM *atomNum,ALLOY *alloy,
						PE energy,PSOPARA *para,char *output);

/***************************************
* @Name:PSO_Start
* @Purpose：运行粒子群算法实例，会在[output]文件夹下生成结果文件和分析文件
* @param：PSOInstance *instance --- 粒子群算法实例
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void PSO_Start(PSOInstance *instance,char *output);

/***************************************
* @Name:updatePbest
* @Purpose：更新个体历史最优
* @param：PSOInstance *instance --- 粒子群算法实例
* @return：void --- 
****************************************/
void updatePbest(PSOInstance *instance);

/***************************************
* @Name:updateGbest
* @Purpose：更新全局最优
* @param：PSOInstance *instance --- 粒子群算法实例
* @return：void --- 
****************************************/
void updateGbest(PSOInstance *instance);

/***************************************
* @Name:updateSwapList
* @Purpose：更新种群中每个个体的交换序
* @param：PSOInstance *instance --- 粒子群算法实例
* @return：void --- 
****************************************/
void updateSwapList(PSOInstance *instance);

/***************************************
* @Name:upDateChorm
* @Purpose：更新种群中每个个体的解序列
* @param：PSOInstance *instance --- 粒子群算法
* @return：void --- 
****************************************/
void upDateChorm(PSOInstance *instance);

/***************************************
* @Name:PSOMutation
* @Purpose：调整种群中每个个体的历史最优解
			调整方式如下：
			交换两个不同类型的原子，如果交换后能量更低，则替换原序列，如果更高，则抛弃
* @param：PSOInstance *instance --- 粒子群算法实例
* @return：void --- 
****************************************/
void PSOMutation(PSOInstance *instance);

/***************************************
* @Name:PSO_PrintMsg
* @Purpose：打印粒子群算法实例的基本信息
* @param：PSOInstance *instance --- 粒子群算法实例
* @return：void --- 
****************************************/
void PSO_PrintMsg(PSOInstance *instance, char *output);

/***************************************
* @Name:PSO_EnergyFile
* @Purpose：输出能量文件，记录每一代有效的最优值
			文件路径为：[output]\\energy.txt
			当输出结束之后，需要手动fclose(fp)
* @param：PSOInstance *instance --- 粒子群算法实例
            FILE **fp --- 文件指针的地址
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void PSO_EnergyFile(PSOInstance *instance, FILE **fp, char *output);

/***************************************
* @Name:PSO_ResultFile
* @Purpose：在在文件夹路径下生成结果文件:
			结果文件 - [output]\\result.txt
			绘图文件 - [output]\\Diamond.txt
* @param：PSOInstance *instance --- 粒子群算法实例
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void PSO_ResultFile(PSOInstance *instance, char *output);

/***************************************
* @Name:PSOIndividual_Init
* @Purpose：初始化个体实例，为其开辟所需的空间
			每一个个体实例都需要初始化
* @param：PSOINDIVIDUAL *one --- 粒子群算法个体实例
            int N --- 总原子数
* @return：void --- 
****************************************/
void PSOIndividual_Init(PSOINDIVIDUAL *one,int N);

/***************************************
* @Name:changePSOIndividual
* @Purpose：将源个体实例中的能量和解序列负责给目标个体实例
			用于将个体实例的值赋值给个体历史最优实例
* @param：PSOINDIVIDUAL *one --- 目标个体实例
            PSOINDIVIDUAL *two --- 源个体实例
            int N --- 原子总个数
* @return：void --- 
****************************************/
void changePSOIndividual(PSOINDIVIDUAL *one,PSOINDIVIDUAL *two,int N);

/***************************************
* @Name:PSOIndividual_Free
* @Purpose：释放个体实例初始化中开辟的空间
			当一个个体实例one运行 PSOIndividual_Init(one) 后，当不使用时候必须执行 PSOIndividual_Free(one)
* @param：PSOINDIVIDUAL *one --- 粒子群算法个体实例
* @return：void --- 
****************************************/
void PSOIndividual_Free(PSOINDIVIDUAL *one);

/***************************************
* @Name:printfPSOIndividual
* @Purpose：打印粒子群个体实例的信息
* @param：PSOINDIVIDUAL *individual --- 粒子群算法个体实例
            int N --- 原子总个数
* @return：void --- 
****************************************/
void printfPSOIndividual(PSOINDIVIDUAL *individual,int N);

/***************************************
* @Name:initSwapSequence
* @Purpose：产生一个长度为size的随机交换序
* @param：int size --- 交换子个数
            int N --- 原子总个数
* @return：PSOSWAP* --- 交换序的首地址
****************************************/
PSOSWAP* initSwapSequence(int size,int N);

/***************************************
* @Name:creatSwap
* @Purpose：产生一个值为element1,element2的交换子
* @param：int element1 --- 交换子值1
            int element2 --- 交换子值2
* @return：PSOSWAP* --- 交换子地址
****************************************/
PSOSWAP* creatSwap(int element1,int element2);

/***************************************
* @Name:swapSubtraction
* @Purpose：解序列的减法
* @param：PSOINDIVIDUAL *best --- 最优个体实例
            PSOINDIVIDUAL *pop --- 个体实例
            int N --- 原子总个数
* @return：PSOSWAP* ---  交换序首地址
****************************************/
PSOSWAP* swapSubtraction(PSOINDIVIDUAL *best, PSOINDIVIDUAL *pop, int N);

/***************************************
* @Name:swapAddition
* @Purpose：交换序的加法
* @param：double w --- 惯性系数
            PSOINDIVIDUAL *pop --- 个体实例 
            PSOSWAP *swapList1 --- 交换序1（历史最优－个体）
            PSOSWAP *swapList2 --- 交换序2（全局最优－个体）
* @return：PSOSWAP* --- 交换序首地址
****************************************/
PSOSWAP* swapAddition(double w, PSOINDIVIDUAL *pop, PSOSWAP *swapList1, PSOSWAP *swapList2);

/***************************************
* @Name:swapMultiplication
* @Purpose：交换序的乘法
* @param：double rate ---　交换子被选择的概率
            PSOSWAP *swapList --- 交换序
            PSOSWAP *newSwapList --- 新交换序
* @return：PSOSWAP* --- 交换序的尾地址
****************************************/
PSOSWAP* swapMultiplication(double rate, PSOSWAP *swapList, PSOSWAP *newSwapList);

/***************************************
* @Name:printfPSOSwap
* @Purpose：打印交换序信息
* @param：PSOSWAP *head --- 交换序
* @return：void --- 
****************************************/
void printfPSOSwap(PSOSWAP *head);


//以下未调试，不可使用
void PSO2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output);
void PSO2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output);
void PSO2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int Popsize,double w,double rate,char *output);
