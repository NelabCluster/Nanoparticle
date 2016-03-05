#pragma once
#include "Base.h"
      
double Rmax;
double Rmin;

ATOMPARA GetAtomPara(ATOM atom);

PEnergy3 GetEnergyFunction(PE type);
PCutEnergy3 GetCutEnergyFunction(PE type);
PEnergy GetEnergyFunction1(PE type);

/***************************************
* @Name:ReadCood
* @Purpose��������״ԭ������ȡ����x,y,z
* @param��char *shape --- ��״ 
            int N --- ��ԭ����
            double *x --- x������ ָ�룬��ǰ���ٿռ�
            double *y --- y������ ָ�룬��ǰ���ٿռ�
            double *z --- z������ ָ�룬��ǰ���ٿռ�
* @return��void --- 
****************************************/
void ReadCood(char *shape,int N,double *x,double *y,double *z);
void ReadCood1(char *shape,int N,COOD *cood);

/***************************************
* @Name:ReadFile
* @Purpose�������ļ���ȡԭ�ӷֲ�note������x,y,z
* @param��char *Input --- �ļ�����
            int *note --- ԭ�ӷֲ� ָ�룬��ǰ���ٿռ�
            double *x --- x������ ָ�룬��ǰ���ٿռ�
            double *y --- y������ ָ�룬��ǰ���ٿռ�
            double *z --- z������ ָ�룬��ǰ���ٿռ�
            int N --- ��ԭ����
* @return��void --- 
****************************************/
void ReadFile(char *Input,int *note,double *x,double *y,double *z,int N);

//��ʼ����
void MixNoteInt3(int *note,int N,int A,int B);
void MixNoteInt(int *note, int N, ATOMNUM *atomNum);

void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//�������(x[],y[],z[],��������R,ԭ�Ӹ���)
void Distance(double *x,double *y,double *z,double *R,int N);
void Distance1(COOD *cood,COODDIS *dis);
//���ش洢��Ե�ַ
char* StoragePath(char *shape,int N,ALLOY *alloy,ATOMNUM *atomNum,char *Output);
//ԭ�����ڲ�
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);

void printResult3(int *note,double *x,double *y,double *z,int N,double bili,char *path);
void printResult(int *note,int N,COOD *cood,char *path);
void printDiamond3(int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,int N,char *path);
void printDiamond(int *note,int N,COOD *cood,ALLOY *alloy,char *path);

void printData(char *Line_Date,char *Line_End,int N);
//����ԭ�������ľ�����������
int* orderCoodFromCore(double *x,double *y,double *z,int N);
int* OrderCoodAlongX(double *x,double *y,double *z,int N);
void L0NoteInt2(int *note,int N,int A,double *x,double *y,double *z);

/***************************************
* @Name:CoodNum3
* @Purpose����λ���ϸ�ԭ�ӵ�ͳ��
* @param��char *input --- �����ļ�
            int N --- ��ԭ����
            char *Output --- �ļ��洢·��
* @return��void --- 
****************************************/
void CoodNum3(char *input,int N,char *Output);


/***************************************
* @Name:pairCorrelation3
* @Purpose������ֲ�����
* @param��char *input --- �����ļ����� 
            int N --- ��ԭ����
            char *Output --- �ļ��洢·��
* @return��void --- 
****************************************/
void pairCorrelation3(char *input,int N,char *Output);

/***************************************
* @Name:ShellCount3
* @Purpose������ԭ��ͳ��
* @param��char *input --- �����ļ�����
            int N --- ��ԭ����
            char *Output --- ����ļ�����
* @return��void --- 
****************************************/
void ShellCount3(char *input,int N,char *Output);

//����JA,JB,�����Ƕ��Ͻ�,�����ӡ����.
void JAJB(char *input,int N);

/***************************************
* @Name:saveMatrix
* @Purpose�������ά�����ڱ���
* @param��double *x --- 
            int N --- 
            char* output --- 
* @return��void --- 
****************************************/
void saveMatrix(double *x, int N, char* output);

/***************************************
* @Name:getLatticeParameter3
* @Purpose���������Ͻ�ľ���
			���������Ͻ�ֻҪ��atom2 = atom3 ����
* @param��ATOM atom1 --- ԭ������1
            ATOM atom2 --- ԭ������2
            ATOM atom3 --- ԭ������3
* @return��double --- ������
****************************************/
double getLatticeParameter3(ATOM atom1,ATOM atom2,ATOM atom3);

/***************************************
* @Name:getLatticeParameter
* @Purpose������ԭ�����ʹ������ĺϽ𾧸�
* @param��ATOM* atoms --- ԭ������
            int atomTypeCount --- ԭ��
* @return��double --- ����
****************************************/
double getLatticeParameter(ALLOY *alloy);

void SD_File(char *input, ATOM *atoms, int atomTypeCount, int N);

void Alloy_Init(ALLOY *alloy,...);
void Alloy_Copy(ALLOY *to, ALLOY *from);
void Alloy_Free(ALLOY *alloy);
void AtomNum_Init(ATOMNUM *atomNum,...);
void AtomNum_Copy(ATOMNUM *to, ATOMNUM *from);
void AtomNum_Free(ATOMNUM *atomNum);

void Cood_Init(COOD *cood,int N);
void Cood_Free(COOD *cood);