#pragma once
/**************************势能原子***********************************/
//	QSC: Ni Cu Rh Pd Ag Ir Pt Au  
//	QSCNewModel: 占时不用
//	AEAM: 占时不用
//	Johnson: Fe Pt Cu 占时不用
//	TBM(Gupta): Pt Pd1 Pd2 Pd3 占时不用
/**************************形状尺寸***********************************/
//	SP:	19 141 459 1061 2171 3781 5895 8685 12479	
//	CU: 63 365 1099 2457 4631 7813 12195
//	RD: 19 93 279 617 1163 1957 3055 4497 6339 8621 11399
//	TH: 35 165 455 969 1771 2925 4495 6545 9139 12341
//	OH: 19 85 231 489 891 1469 2255 3281 4579 6181 8119 10425
//	TO: 63 309 711 1289 2075 3101 4399 6001 7939 10245
//	DH: 104 608 1832 4096 7720 13024
//	IH: 308 2056 6524 14992
//	THH210: 63 443 1417 3285 6323 10831 25385 49251
//	THH310: 1249 9553 31825 74977
#include "Base.h"		//基本数据的定义
#include "Tool.h"		//通用函数的定义
#include "Algorithm.h"	//算法的定义
#include "Energy.h"		//势能的定义

main()
{ 	

	/*
		以遗传算法为例，说明如何使用：
		GAPARA para：定义算法参数实例，该参数包含的内容可在GA.H中查看，下面一一列出，
					 	int popSize：				//种群个数
						double pc：					//交叉概率
						double pm;					//变异概率
						int maxGenerations;			//最大迭代次数 (当迭代次数达到最大迭代次数 
						int convergenceGenerations;	//收敛代数		或者 不变的代数超过收敛代数时候， 算法终止)
						char needOrderCood;			//是否需要排序；1：需要，0：不需要
						char needAdjustment;		//是否需要调整算子；1：需要，0：不需要
						double adjustmentRate;		//调整概率，当needAdjustment是1时有效
		刚定义的参数实例未赋值，可以通过一下函数赋默认值：
		GAPara_Init( &para ); 默认值在GA.H中可以查看，这里列出：
								int popSize = 100：			
								double pc = 0.8：				
								double pm = 0.0015;				
								int maxGenerations = 5000;			
								int convergenceGenerations = 1000;	
								char needOrderCood = 1;			
								char needAdjustment = 1;		
								double adjustmentRate = 1.0;
		也可以手动修改，如设置最大迭代次数为3000代：
		para.maxGeneration = 3000;
		设置好参数后，就可以带入二合金或是三合金的函数中进行运算：
		GA2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"GA2"); 
		GA3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"GA3");
		二合金每个参数的意义如下：	"CU":构型是立方体
									63	:原子总数是63
									30	:Pt的原子数是30
									Pt,Pd:原子类型为Pt,Pd
									QSC	:势能为QSC
									&para:遗传算法参数
									"GA2":结果存在的文件夹名称GA2
		三合金每个参数的意义如下：	"CU":构型是立方体
									63	:原子总数是63
									10,20	:Pt的原子数是10,Pd原子个数是20
									Pt,Pd,Au:原子类型为Pt,Pd,Au
									QSC	:势能为QSC
									&para:遗传算法参数
									"GA3":结果存在的文件夹名称GA3

		1、算法的优化是自动开启的，可以通过将参数中 needOrderCood 和 needAdjustment 设置为0 取消优化；
		   当参数进行修改后，记得修改存储结果的文件夹名称，防止文件的覆盖。
		2、当遇到算法长时间无法收敛的情况，可以通过将参数中的 maxGenerations 或是 convergenceGenerations 设置较小值，降低收敛条件。
		   需要注意的是如果是不同算法对于同一个对象的对比，最好设置相同的收敛条件。
	*/
//	遗传算法
//	GAPARA para;
//	GAPara_Init(&para);
//	GA2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"GA2");
//	GA3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"GA3");

//	蒙特卡罗法
//	MCPARA para;
//	MCPara_Init(&para);
//	MC2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"MC2");
//	MC3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"MC3");

//	粒子群算法
//	PSOPARA para;
//	PSOPara_Init( &para );
//	PSO2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"PSO2");
//	PSO3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"PSO3");
	
//	鱼群算法未完成，不可使用
//	FISHPARA para;
//	FISHPara_Init( &para );
//	FISH2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"FISH2");
//	FISH3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"PSO3");

}