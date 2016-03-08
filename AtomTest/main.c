#pragma once
/**************************����ԭ��***********************************/
//	QSC: Ni Cu Rh Pd Ag Ir Pt Au  
//	QSCNewModel: ռʱ����
//	AEAM: ռʱ����
//	Johnson: Fe Pt Cu ռʱ����
//	TBM(Gupta): Pt Pd1 Pd2 Pd3 ռʱ����
/**************************��״�ߴ�***********************************/
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
#include "Base.h"		//�������ݵĶ���
#include "Tool.h"		//ͨ�ú����Ķ���
#include "Algorithm.h"	//�㷨�Ķ���
#include "Energy.h"		//���ܵĶ���

main()
{ 	

	/*
		���Ŵ��㷨Ϊ����˵�����ʹ�ã�
		GAPARA para�������㷨����ʵ�����ò������������ݿ���GA.H�в鿴������һһ�г���
					 	int popSize��				//��Ⱥ����
						double pc��					//�������
						double pm;					//�������
						int maxGenerations;			//���������� (�����������ﵽ���������� 
						int convergenceGenerations;	//��������		���� ����Ĵ���������������ʱ�� �㷨��ֹ)
						char needOrderCood;			//�Ƿ���Ҫ����1����Ҫ��0������Ҫ
						char needAdjustment;		//�Ƿ���Ҫ�������ӣ�1����Ҫ��0������Ҫ
						double adjustmentRate;		//�������ʣ���needAdjustment��1ʱ��Ч
		�ն���Ĳ���ʵ��δ��ֵ������ͨ��һ�º�����Ĭ��ֵ��
		GAPara_Init( &para ); Ĭ��ֵ��GA.H�п��Բ鿴�������г���
								int popSize = 100��			
								double pc = 0.8��				
								double pm = 0.0015;				
								int maxGenerations = 5000;			
								int convergenceGenerations = 1000;	
								char needOrderCood = 1;			
								char needAdjustment = 1;		
								double adjustmentRate = 1.0;
		Ҳ�����ֶ��޸ģ�����������������Ϊ3000����
		para.maxGeneration = 3000;
		���úò����󣬾Ϳ��Դ�����Ͻ�������Ͻ�ĺ����н������㣺
		GA2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"GA2"); 
		GA3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"GA3");
		���Ͻ�ÿ���������������£�	"CU":������������
									63	:ԭ��������63
									30	:Pt��ԭ������30
									Pt,Pd:ԭ������ΪPt,Pd
									QSC	:����ΪQSC
									&para:�Ŵ��㷨����
									"GA2":������ڵ��ļ�������GA2
		���Ͻ�ÿ���������������£�	"CU":������������
									63	:ԭ��������63
									10,20	:Pt��ԭ������10,Pdԭ�Ӹ�����20
									Pt,Pd,Au:ԭ������ΪPt,Pd,Au
									QSC	:����ΪQSC
									&para:�Ŵ��㷨����
									"GA3":������ڵ��ļ�������GA3

		1���㷨���Ż����Զ������ģ�����ͨ���������� needOrderCood �� needAdjustment ����Ϊ0 ȡ���Ż���
		   �����������޸ĺ󣬼ǵ��޸Ĵ洢������ļ������ƣ���ֹ�ļ��ĸ��ǡ�
		2���������㷨��ʱ���޷����������������ͨ���������е� maxGenerations ���� convergenceGenerations ���ý�Сֵ����������������
		   ��Ҫע���������ǲ�ͬ�㷨����ͬһ������ĶԱȣ����������ͬ������������
	*/
//	�Ŵ��㷨
//	GAPARA para;
//	GAPara_Init(&para);
//	GA2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"GA2");
//	GA3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"GA3");

//	���ؿ��޷�
//	MCPARA para;
//	MCPara_Init(&para);
//	MC2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"MC2");
//	MC3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"MC3");

//	����Ⱥ�㷨
//	PSOPARA para;
//	PSOPara_Init( &para );
//	PSO2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"PSO2");
//	PSO3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"PSO3");
	
//	��Ⱥ�㷨δ��ɣ�����ʹ��
//	FISHPARA para;
//	FISHPara_Init( &para );
//	FISH2_InitWithMixing("CU",63,30,Pt,Pd,QSC,&para,"FISH2");
//	FISH3_InitWithMixing("CU",63,10,20,Pt,Pd,Au,QSC,&para,"PSO3");

}