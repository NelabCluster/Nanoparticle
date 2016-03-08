#pragma once
#include "Base.h"
#include "Energy.h"
#include "Tool.h"

/***************************************
* @Name: 蒙特卡罗法算法参数
****************************************/
typedef struct
{
	int RSTEP;			//最大迭代次数
	int temperature;	//温度，保留，未使用
} MCPARA;

static const MCPARA _defMCPara = {
	100000,0
};

/***************************************
* @Name: 蒙特卡罗法算法个体
****************************************/
typedef struct
{
	int *note;			//原子序列
	double E;			//能量
} MCINDIVIDUAL;

/***************************************
* @Name: 蒙特卡罗法算法实例
****************************************/
typedef struct
{
	char *shape;		//构型名称	
	int N;				//原子总数
	PE energyType;		//势能类型
	ALLOY alloy;		//合金类型
	ATOMNUM atomNum;	//各类金属的原子个数
	COOD cood;			//坐标值
	COODDIS dis;		//原子间距离
	MCPARA para;		//算法参数
	MCINDIVIDUAL one;	//个体实例	
	int clocks;			//迭代次数
	int IJKL;			//找到最优值的次数
} MCInstance;


/***************************************
* @Name:MCPara_Init
* @Purpose：给para赋值默认的值
			默认值为静态变量 _defMCPara
* @param：MCPARA *para --- 蒙特卡罗算法参数
* @return：void --- 
****************************************/
void MCPara_Init(MCPARA *para);

/***************************************
* @Name:MC2_InitWithMixing
* @Purpose：蒙特卡罗二合金初始结构为随机分布
* @param：char *shape --- 形状
            int N --- 总原子数
            int A --- 原子A个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
			PE energy --- 势能
            MCPARA *para --- 蒙特卡罗法算法参数
            char *Output --- 存储文件名
* @return：void --- 
****************************************/
void MC2_InitWithMixing(char *shape,int N,int A,
						ATOM atomA,ATOM atomB,
						PE energy,MCPARA *para,char *output);

/***************************************
* @Name:MC3_InitWithMixing
* @Purpose：蒙特卡罗算法 三合金 初始结构为随机分布
* @param：char *shape --- 形状
            int N --- 总原子数
            int A --- 原子A个数
            int B --- 原子B个数
            ATOM atomA --- 原子A类型
            ATOM atomB --- 原子B类型
            ATOM atomC --- 原子C类型
			PE energy --- 势能
            MCPARA *para --- 蒙特卡罗法算法参数
            char *Output --- 存储文件夹名称
* @return：void --- 
****************************************/
void MC3_InitWithMixing(char *shape,int N,
						int A,int B,
						ATOM atomA,ATOM atomB,ATOM atomC,
						PE energy,MCPARA *para,char *output);

/***************************************
* @Name:MC_InitWithMixing
* @Purpose：蒙特卡罗算法 多合金 初始结构随机排布
* @param：char *shape --- 构型名称
            int N --- 原子总个数
            ATOMNUM *atomNum --- 各类金属原子数量 
            ALLOY *alloy --- 合金类型
            PE energy --- 能量类型
            MCPARA *para --- 蒙特卡罗算法参数
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void MC_InitWithMixing(char *shape, int N, 
					   ATOMNUM *atomNum, ALLOY *alloy, 
					   PE energy, MCPARA *para, char *output);

/***************************************
* @Name:MC_Start
* @Purpose：运行蒙特卡罗算法实例，会在[output]文件夹下生成结果文件和分析文件
* @param：MCInstance *instance --- 蒙特卡罗算法实例
            char *output --- 存储文件夹名称
* @return：void --- 
****************************************/
void MC_Start(MCInstance *instance, char *output);


/***************************************
* @Name:MC_PrintMsg
* @Purpose：打印蒙特卡罗算法实例的信息
* @param：MCInstance *instance --- 蒙特卡罗算法实例
* @return：void --- 
****************************************/
void MC_PrintMsg( MCInstance *instance, char *output );

/***************************************
* @Name:MC_EnergyFile
* @Purpose：输出能量文件，记录每一代有效的最优值
			文件路径为：[output]\\energy.txt
			当输出结束之后，需要手动fclose(fp)
* @param：MCInstance *instance --- 蒙特卡罗算法实例
            FILE **fp --- 文件指针的地址
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void MC_EnergyFile(MCInstance *instance, FILE **fp, char *output);

/***************************************
* @Name:MC_ResultFile
* @Purpose：在在文件夹路径下生成结果文件:
			结果文件 - [output]\\result.txt
			绘图文件 - [output]\\Diamond.txt
* @param：MCInstance *instance --- 蒙特卡罗算法实例
            char *output --- 文件夹路径
* @return：void --- 
****************************************/
void MC_ResultFile(MCInstance *instance, char *output);

/***************************************
* @Name:MCIndividual_Init
* @Purpose：初始化个体实例，为其开辟所需的空间
			每一个个体实例都需要初始化
* @param：MCINDIVIDUAL *one --- 蒙特卡罗算法个体实例
            int N --- 总原子数
* @return：void --- 
****************************************/
void MCIndividual_Init(MCINDIVIDUAL *one,int N);

/***************************************
* @Name:MCIndividual_Free
* @Purpose：释放个体实例初始化中开辟的空间
			当一个个体实例one运行 GAIndividual_Init(one) 后，
			当不使用时候必须执行 GAIndividual_Free(one)
* @param：MCINDIVIDUAL *one --- 蒙特卡罗算法个体实例
* @return：void --- 
****************************************/
void MCIndividual_Free(MCINDIVIDUAL *one);



//以下未调试,不可使用
//初始相位分离
void MC2_InitWithPhase(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
//初始L0结构，可能不完整
void MC2_InitWithL0(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);

//初始核壳，壳层的层数为surfaceLayer 0:最外层
void MC2_InitWithFromCore(char *shape,int N,int A,ATOM atomA,ATOM atomB,PE energy,int RSTEP,char *Output);
void MC3_InitWithCoreSurface(char *shape,int N,int surfaceLayer,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
void MC3_InitWithFromCore(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,PE energy,int RSTEP,char *Output);
