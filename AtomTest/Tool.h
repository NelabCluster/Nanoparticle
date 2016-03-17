#pragma once
#include "Base.h"
#include "Energy.h"

double Rmax;
double Rmin;

ATOMPARA GetAtomPara(ATOM atom);

void Energy_Init( PE type, ALLOY *alloy );
void Energy_Free( PE type );
PEnergy GetEnergyFunction(PE type);
PCutEnergy GetCutEnergyFunction(PE type);

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
void ReadFile1( char *input, int *note, COOD *cood, int N );

//��ʼ����
void MixNoteInt(int *note, int N, ATOMNUM *atomNum);

void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//�������(x[],y[],z[],��������R,ԭ�Ӹ���)
void Distance(double *x,double *y,double *z,double *R,int N);
void Distance1(COOD *cood,COODDIS *dis);

/***************************************
* @Name:StoragePath
* @Purpose�����ݹ��ͳߴ�Ͻ����ͺ͸������ļ�����������Ӧ��·��
			·����[Output]\\��������\\ԭ������\\�Ͻ�����\\ԭ�ӱ���
			�ڸ�·���´洢����ļ��ͷ����ļ�
* @param��char *shape --- ��������
            int N --- ԭ������
            ALLOY *alloy --- �Ͻ�����
            ATOMNUM *atomNum --- �����������
            char *Output --- �ļ�������
* @return��char* --- 
****************************************/
char* StoragePath(char *shape,int N,ALLOY *alloy,ATOMNUM *atomNum,char *Output);

//ԭ�����ڲ�
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);


/***************************************
* @Name:printResult
* @Purpose�����ļ�·�������������ļ�
			note����Ҫ���ٹ��ռ�
			cood: ��Ҫ��ʼ����
			path��������txt�ļ���һ��Ϊ ...\\result.txt
* @param��int *note --- ����
            int N --- ԭ������
            COOD *cood --- ����
            char *path --- �ļ���·��
* @return��void --- 
****************************************/
void printResult( int *note, int N, COOD *cood, char *path );

/***************************************
* @Name:printDiamond
* @Purpose�����ļ�·�����������ͼ�ļ�
			note����Ҫ���ٹ��ռ�
			cood����Ҫ��ʼ����
			alloy����Ҫ��ʼ����
			path��������txt�ļ���һ��Ϊ ...\\Diamond.txt
* @param��int *note --- ����
            int N --- ԭ������
            COOD *cood --- ����
            ALLOY *alloy --- �Ͻ�
            char *path --- �ļ���·��
* @return��void --- 
****************************************/
void printDiamond(int *note,int N,COOD *cood,ALLOY *alloy,char *path);

void printData(char *Line_Date, int N);

//����ԭ�������ľ�����������
int* orderCoodFromCore(double *x,double *y,double *z,int N);
int* OrderCoodAlongX(double *x,double *y,double *z,int N);
void L0NoteInt2(int *note,int N,int A,double *x,double *y,double *z);

/***************************************
* @Name:CoordinateNumber
* @Purpose����λ���ϸ�ԭ�ӵ�ͳ��
* @param��char *input --- �����ļ�
            int N --- ��ԭ����
            char *output --- �ļ��洢·��
* @return��void --- 
****************************************/
void CoordinateNumber( char *input, int N, char *output );

/***************************************
* @Name:coordinationNumberR
* @Purpose�����ݾ��뷵�ظ���ԭ�ӵ���λ��
			������ֵ��ʹ��ʱ����Ҫ�ֶ��ͷ�
* @param��COODDIS *dis --- ����
* @return��int* --- ��λ������
****************************************/
int* coordinationNumberR( COODDIS *dis );

/***************************************
* @Name:pairCorrelation
* @Purpose������ֲ�����
* @param��char *input --- �����ļ����� 
            int N --- ��ԭ����
            char *Output --- �ļ��洢·��
* @return��void --- 
****************************************/
void pairCorrelation(char *input,int N,char *Output);

/***************************************
* @Name:ShellCount
* @Purpose������ԭ��ͳ��
* @param��char *input --- �����ļ�����
            int N --- ��ԭ����
            char *Output --- ����ļ�����
* @return��void --- 
****************************************/
void ShellCount( char *input, int N, char *output );

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