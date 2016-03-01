#pragma once
#include "Base.h"
      
double Rmax;
double Rmin;

ATOMPARA GetAtomPara(ATOM atom);

PEnergy3 GetEnergyFunction(PE type);
PCutEnergy3 GetCutEnergyFunction(PE type);


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
void MixNoteInt(int *note, int N, ATOMNUM atomNum);

void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//�������(x[],y[],z[],��������R,ԭ�Ӹ���)
void Distance(double *x,double *y,double *z,double *R,int N);
void Distance1(COOD cood,COODDIS *dis);
//���ش洢��Ե�ַ
char* StoragePath(char *shape,int N,int A,int B,ATOM atomA,ATOM atomB,ATOM atomC,char *Output);
//ԭ�����ڲ�
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);

void printResult3(int *note,double *x,double *y,double *z,int N,double bili,char *path);
void printResult(int *note,double *x,double *y,double *z,int N,double bili,char *path);
void printDiamond3(int *note,double *x,double *y,double *z,ATOM atomA,ATOM atomB,ATOM atomC,int N,char *path);
void printDiamond(int *note,double *x,double *y,double *z,ALLOY alloy,int N,char *path);

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
double getLatticeParameter(ATOM* atoms, int atomTypeCount);

void SD_File(char *input, ATOM *atoms, int atomTypeCount, int N);