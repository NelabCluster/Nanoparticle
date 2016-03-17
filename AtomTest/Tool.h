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
* @Purpose：根据形状原子数读取坐标x,y,z
* @param：char *shape --- 形状 
            int N --- 总原子数
            double *x --- x轴坐标 指针，提前开辟空间
            double *y --- y轴坐标 指针，提前开辟空间
            double *z --- z轴坐标 指针，提前开辟空间
* @return：void --- 
****************************************/
void ReadCood(char *shape,int N,double *x,double *y,double *z);
void ReadCood1(char *shape,int N,COOD *cood);

/***************************************
* @Name:ReadFile
* @Purpose：根据文件读取原子分布note和坐标x,y,z
* @param：char *Input --- 文件名称
            int *note --- 原子分布 指针，提前开辟空间
            double *x --- x轴坐标 指针，提前开辟空间
            double *y --- y轴坐标 指针，提前开辟空间
            double *z --- z轴坐标 指针，提前开辟空间
            int N --- 总原子数
* @return：void --- 
****************************************/
void ReadFile(char *Input,int *note,double *x,double *y,double *z,int N);
void ReadFile1( char *input, int *note, COOD *cood, int N );

//初始排列
void MixNoteInt(int *note, int N, ATOMNUM *atomNum);

void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B);
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z);
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z);
void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z);
//计算距离(x[],y[],z[],距离数组R,原子个数)
void Distance(double *x,double *y,double *z,double *R,int N);
void Distance1(COOD *cood,COODDIS *dis);

/***************************************
* @Name:StoragePath
* @Purpose：根据构型尺寸合金类型和个数在文件夹下生成相应的路径
			路径：[Output]\\构型名称\\原子总数\\合金类型\\原子比例
			在该路径下存储结果文件和分析文件
* @param：char *shape --- 构型名称
            int N --- 原子总数
            ALLOY *alloy --- 合金类型
            ATOMNUM *atomNum --- 各类金属个数
            char *Output --- 文件夹名称
* @return：char* --- 
****************************************/
char* StoragePath(char *shape,int N,ALLOY *alloy,ATOMNUM *atomNum,char *Output);

//原子所在层
int* Shell_Shape(char *shape,int N);
int* Shell_Cood(double *x,double *y,double *z,int N);


/***************************************
* @Name:printResult
* @Purpose：在文件路径名中输出结果文件
			note：需要开辟过空间
			cood: 需要初始化过
			path：必须是txt文件，一般为 ...\\result.txt
* @param：int *note --- 序列
            int N --- 原子总数
            COOD *cood --- 坐标
            char *path --- 文件名路径
* @return：void --- 
****************************************/
void printResult( int *note, int N, COOD *cood, char *path );

/***************************************
* @Name:printDiamond
* @Purpose：在文件路径名中输出绘图文件
			note：需要开辟过空间
			cood：需要初始化过
			alloy：需要初始化过
			path：必须是txt文件，一般为 ...\\Diamond.txt
* @param：int *note --- 序列
            int N --- 原子总数
            COOD *cood --- 坐标
            ALLOY *alloy --- 合金
            char *path --- 文件名路径
* @return：void --- 
****************************************/
void printDiamond(int *note,int N,COOD *cood,ALLOY *alloy,char *path);

void printData(char *Line_Date, int N);

//根据原子里质心距离排序坐标
int* orderCoodFromCore(double *x,double *y,double *z,int N);
int* OrderCoodAlongX(double *x,double *y,double *z,int N);
void L0NoteInt2(int *note,int N,int A,double *x,double *y,double *z);

/***************************************
* @Name:CoordinateNumber
* @Purpose：配位数上各原子的统计
* @param：char *input --- 输入文件
            int N --- 总原子数
            char *output --- 文件存储路径
* @return：void --- 
****************************************/
void CoordinateNumber( char *input, int N, char *output );

/***************************************
* @Name:coordinationNumberR
* @Purpose：根据距离返回各个原子的配位数
			当返回值不使用时，需要手动释放
* @param：COODDIS *dis --- 距离
* @return：int* --- 配位数序列
****************************************/
int* coordinationNumberR( COODDIS *dis );

/***************************************
* @Name:pairCorrelation
* @Purpose：径向分布函数
* @param：char *input --- 输入文件名称 
            int N --- 总原子数
            char *Output --- 文件存储路径
* @return：void --- 
****************************************/
void pairCorrelation(char *input,int N,char *Output);

/***************************************
* @Name:ShellCount
* @Purpose：各层原子统计
* @param：char *input --- 输入文件名称
            int N --- 总原子数
            char *Output --- 输出文件名称
* @return：void --- 
****************************************/
void ShellCount( char *input, int N, char *output );

//计算JA,JB,输入是二合金,结果打印出来.
void JAJB(char *input,int N);

/***************************************
* @Name:saveMatrix
* @Purpose：保存二维数组于本地
* @param：double *x --- 
            int N --- 
            char* output --- 
* @return：void --- 
****************************************/
void saveMatrix(double *x, int N, char* output);

/***************************************
* @Name:getLatticeParameter
* @Purpose：计算原子类型大于三的合金晶格
* @param：ATOM* atoms --- 原子类型
            int atomTypeCount --- 原子
* @return：double --- 晶格
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