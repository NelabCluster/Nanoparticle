#pragma once
/**************************����ԭ��***********************************/
//	QSC: Ni Cu Rh Pd Ag Ir Pt Au  
//	QSCNewModel: ռʱ����
//	AEAM:
//	Johnson: Fe Pt Cu
//	TBM(Gupta): Pt Pd1 Pd2 Pd3
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

	//�㷨�����ؿ����㷨��
	//��ʼ�ṹ�����
	//��״��������
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Au):180
	//�Ͻ�����: Au-Ag;
	//���ܣ�QSC
	//������10000
	//�洢���ļ������ƣ���MCS_MIX���������Լ��趨
//	MCS2_InitWithMixing("CU",365,180,Au,Ag,QSC,100000,"MCS_MIX");

	/****������ʼ�ṹ*****/
	//��ʼ�ṹ����λ����
	//�洢���ļ������ƣ���MCS_PHASE��
//	MCS2_InitWithPhase("CU",365,180,Au,Ag,QSC,10000,"MCS_PHASE");
	//��ʼ�ṹ���˿ǽṹ����һ��ԭ�����ڲ�
	//�洢���ļ������ƣ���MCS_CORE��
//	MCS2_InitWithFromCore("CU",365,180,Au,Ag,QSC,10000,"MCS_CORE");
	//��ʼ�ṹ: �����ʼ�ṹ
	//�洢���ļ������ƣ���MCS_L0��
//	MCS2_InitWithL0("CU",365,180,Au,Ag,QSC,10000,"MCS_L0");


	//�㷨���Ŵ��㷨��
	//��ʼ�ṹ�����
	//��״�������� CU
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Au):180
	//�Ͻ�����: Au-Ag
	//���ܣ�QSC
	//��Ⱥ��С��100
	//������ʣ�0.8
	//�����ʣ�0.0015
	//�������ʣ�1.0
	//�洢���ļ������ƣ�GA2_MIX
//	GA2_InitWithMixing("CU",365,180,Au,Ag,QSC,100,0.8,0.0015,1.0,"GA2_MIX");

	//�㷨���Ŵ��㷨��
	//��ʼ�ṹ�����
	//��״�������� CU
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Pt):100
	//�ڶ���ԭ�ӵĸ���(�������Pd):200
	//�Ͻ�����: Pt-Pd-Au
	//���ܣ�QSC
	//��Ⱥ��С��100
	//������ʣ�0.8
	//�����ʣ�0.0015
	//�������ʣ�1.0
	//�洢���ļ������ƣ�GA3_MIX
//	GA3_InitWithMixing("CU",365,100,200,Pt,Pd,Au,QSC,100,0.8,0.0015,1.0,"GA3_MIX");

	//�㷨������Ⱥ�㷨��
	//��ʼ�ṹ�����
	//��״��������
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Pt):180
	//�Ͻ�����: Au-Ag
	//���ܣ�QSC
	//��Ⱥ��С��100
	//����ϵ����0.25
	//�����ʣ�1.0
//	PSO2_InitWithMixing("CU",365,180,Au,Ag,QSC,100,0.25,1.0,"PSO2_MIX_YES",PSOYesImprove);//�Ľ���PSO �洢���ļ������ƣ�PSO_PARA1_YES
//	PSO2_InitWithMixing("CU",365,180,Au,Ag,QSC,100,0.25,1.0,"PSO2_MIX_NO",PSONoImprove);//δ�Ľ�PSO //�洢���ļ������ƣ�PSO_PARA1_NO
	

	//�㷨������Ⱥ�㷨��
	//��ʼ�ṹ�����
	//��״�������� CU
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Pt):100
	//�ڶ���ԭ�ӵĸ���(�������Pd):200
	//�Ͻ�����: Pt-Pd-Au
	//���ܣ�QSC
	//���ܣ�QSC
	//��Ⱥ��С��100
	//����ϵ����0.25
	//�����ʣ�1.0
//	PSO3_InitWithMixing("CU",365,100,200,Pt,Pd,Au,QSC,100,0.25,1.0,"PSO3_MIX");//�Ľ���PSO �洢���ļ������ƣ�PSO_PARA1_YES


}