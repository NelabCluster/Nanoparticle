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

void test()
{
	double *x,*y,*z,*R;
	int *note;
	int N;
	double E1,E2,E3;
	ATOM atoms[3] = {Pt,Pd,Rh};

	N = 1417;
	note = calloc(N,sizeof(int));
	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	R = calloc(N*N,sizeof(double));
	ReadFile("result.txt",note,x,y,z,N);
	Distance(x,y,z,R,N);
	SetEnergyPow(Pt,Pd,Rh);
	E1 = QSCEnergy3(note,R,Pt,Pd,Rh,N);
	E2 = QSCEnergy(note,R,atoms,3,N);
	printf("%f\n%f\n%f\n",E1,E2,E3);
	system("pause");
}

void test1()
{
	ATOM atoms[3] = {Pt,Pd,Rh};
	int N = 1417;

	SD_File("result.txt",atoms,3,N);
}

main()
{ 
	test1();

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
	//��״��������
	//��ԭ����: 365
	//��һ��ԭ�ӵĸ���(�������Au):180
	//�Ͻ�����: Au-Ag
	//���ܣ�QSC
	//��Ⱥ��С��100
	//������ʣ�0.8
	//�����ʣ�0.0015
	//�������ʣ�1.0
	//�洢���ļ������ƣ���GA_MIX��
//	GA2_InitWithMixing("CU",365,180,Au,Ag,QSC,100,0.8,0.0015,1.0,"GA_MIX");

	//�㷨������Ⱥ�㷨��
	//��ʼ�ṹ�����
	//��״��24����
	//��ԭ����: 3285
	//��һ��ԭ�ӵĸ���(�������Pt):1642
	//�Ͻ�����: Pt-Pd
	//���ܣ�TBM
	//��Ⱥ��С��100
	//����ϵ����0.25
	//�����ʣ�1.0
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd1,TBM,100,0.25,1.0,"PSO_PARA1_YES",PSOYesImprove);//�Ľ���PSO����һ�����
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd2,TBM,100,0.25,1.0,"PSO_PARA2_YES",PSOYesImprove);//�Ľ���PSO���ڶ������
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd3,TBM,100,0.25,1.0,"PSO_PARA3_YES",PSOYesImprove);//�Ľ���PSO�����������

//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd1,TBM,100,0.25,1.0,"PSO_PARA1_NO",PSONoImprove);//�Ľ�ǰPSO����һ�����
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd2,TBM,100,0.25,1.0,"PSO_PARA2_NO",PSONoImprove);//�Ľ�ǰPSO���ڶ������
//	PSO2_InitWithMixing("THH210",3285,1642,Pt,Pd3,TBM,100,0.25,1.0,"PSO_PARA3_NO",PSONoImprove);//�Ľ�ǰPSO�����������

}