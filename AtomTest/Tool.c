#include "Tool.h"
#include "Energy.h"

ATOMPARA GetAtomPara(ATOM atom)
{
	ATOMPARA zero = {0.0,""};
	switch (atom)
	{
	case Ni:
		return NiPara;
	case Cu:
		return CuPara;
	case Rh:
		return RhPara;
	case Pd:
		return PdPara;
	case Ag:
		return AgPara;
	case Ir:
		return IrPara;
	case Pt:
		return PtPara;
	case Au:
		return AuPara;
	case Fe:
		return FePara;
	case Pd1:
		return Pd1Para;
	case Pd2:
		return Pd2Para;
	case Pd3:
		return Pd3Para;
	default:
		break;
	}
	
	return zero;
}

PEnergy3 GetEnergyFunction(PE type)
{
	switch (type)
	{
	case QSC:
		return QSCEnergy3;
	case QSCNewModel:
		return NULL;
	case AEAM:
		return AEAMEnergy;
	case Johnson:
		return JohnsonEnergy;
	case TBM:
		return TBMEnergy;
		break;
	}
	return NULL;
}

PCutEnergy3 GetCutEnergyFunction(PE type)
{
	switch (type)
	{
	case QSC:
		return QSCCutEnergy3;
	case QSCNewModel:
		return NULL;
	case AEAM:
		return AEAMCutEnergy;
	case Johnson:
		return JohnsonCutEnergy;
	case TBM:
		return TBMCutEnergy;
		break;
	}
	return NULL;
}

void Energy_Init( PE type, ALLOY *alloy )
{
	switch( type )
	{
	case QSC:
		{
			QSC_Init( alloy );
			break;
		}
	}
}
void Energy_Free( PE type )
{
	switch( type )
	{
	case QSC:
		{
			QSC_Free();
			break;
		}
	}
}

PEnergy GetEnergyFunction1(PE type)
{
	switch (type)
	{
	case QSC:
		return QSCEnergy;
		break;
	}
	return NULL;
}

PCutEnergy GetCutEnergyFunction1(PE type)
{
	switch (type)
	{
	case QSC:
		return QSCCutEnergy;
		break;
	}
	return NULL;
}

//读取坐标
void ReadCood(char *shape,int N,double *x,double *y,double *z){
	int i;
	char cood[30] ="cood\\";
	char num[100];
	FILE *fp;
	strcat(cood,shape);
	strcat(cood,"\\");
	strcat(cood,shape);
	sprintf(num,"_%d.txt",N);
	strcat(cood,num);
	fp = fopen(cood,"r");
	fgets(num,99,fp);
	for(i=0;i<N;i++){
		fscanf(fp,"%lf %lf %lf",&x[i],&y[i],&z[i]);
	}
	fclose(fp);
}

void ReadCood1(char *shape,int N,COOD *_cood){
	int i;
	char cood[30] ="cood\\";
	char num[100];
	FILE *fp;
	strcat(cood,shape);
	strcat(cood,"\\");
	strcat(cood,shape);
	sprintf(num,"_%d.txt",N);
	strcat(cood,num);
	fp = fopen(cood,"r");
	fgets(num,99,fp);
	for(i=0;i<N;i++){
		fscanf(fp,"%lf %lf %lf",_cood->x+i,_cood->y+i,_cood->z+i);
	}
	_cood->N = N;
	fclose(fp);
}

void ReadFile(char *Input,int *note,double *x,double *y,double *z,int N){
	int i;
	double bili;
	FILE *fp;

	fp = fopen(Input,"r");
	fscanf(fp,"%d %lf",&N,&bili);
	for(i=0;i<N;i++){
		fscanf(fp,"%d %lf %lf %lf",&note[i],&x[i],&y[i],&z[i]);
	}
	fclose(fp);	
}

void ReadFile1( char *input, int *note, COOD *cood, int N )
{
	int i;
	FILE *fp;

	fp = fopen( input, "r" );
	fscanf( fp, "%d", &N );
	for(i=0;i<N;i++){
		fscanf( fp,"%d %lf %lf %lf", &note[i], &cood->x[i], &cood->y[i], &cood->z[i]);
	}
	fclose(fp);	
}

void MixNoteInt(int *note,int N,ATOMNUM* atomNum){
	int i;
	int point,temp,index;
	int *sumNumberBeforeAtom;

	sumNumberBeforeAtom = malloc(atomNum->atomTypeCount * sizeof(int));
	sumNumberBeforeAtom[0] = atomNum->numberOfAtom[0];

	for(i=1;i<atomNum->atomTypeCount;i++)
		sumNumberBeforeAtom[i] = sumNumberBeforeAtom[i-1] + atomNum->numberOfAtom[i];

	for(i=0;i<N;i++)
		note[i] = i;
	for(i=0;i<N;i++){
		point = rand()%(N-i);
		temp = note[i];
		note[i] = note[point+i];
		note[point+i] = temp;

		index = 0;
		while( note[i] >= sumNumberBeforeAtom[index] )
			index ++;
		note[i] = index;
	}

	free(sumNumberBeforeAtom);
}

//完全核壳的初始结构
void CoreSurfaceNote3(int *note,int N,double *x,double *y,double *z,int surfaceLayer,int B)
{	
	int i,j;
	int *shell;
	int coreNumber = 0,surfaceNumber = 0;
	int *surfaceIndex;
	int *surfaceCood;
	int point;
	int temp;

	shell = Shell_Cood(x,y,z,N);

	for(i = 0;i < N;i++)
	{
		if(shell[i] > surfaceLayer) 
		{
			note[i] = 0;
			coreNumber ++;
		}
	}

	surfaceNumber = N - coreNumber;
	surfaceIndex = calloc(surfaceNumber,sizeof(int));
	surfaceCood = calloc(surfaceNumber,sizeof(int));

	j = 0;
	for(i=0;i<N;i++)
	{
		if(shell[i] <= surfaceLayer)
		{
			surfaceIndex[j] = i;
			j ++;
		}
	}

	for(i=0;i<surfaceNumber;i++)
		surfaceCood[i] = i;
	for(i=0;i<surfaceNumber;i++){
		point = rand() % (surfaceNumber-i);
		temp = surfaceCood[i];
		surfaceCood[i] = surfaceCood[point+i];
		surfaceCood[point+i] = temp;
		if(surfaceCood[i] < B)
			surfaceCood[i] = 1;
		else 
			surfaceCood[i] = 2;
	}

	for(i = 0;i < surfaceNumber;i ++)
		note[surfaceIndex[i]] = surfaceCood[i];
}


//从里到外
void FromCoreNoteInt3(int *note,int N,int A,int B,double *x,double *y,double *z){
	int i,j;
	int point;
	int temp;
	int coreNumber = 0,surfaceNumber = 0;
	int *surfaceIndex;
	int *surfaceCood;
	int *index;

	index = orderCoodFromCore(x,y,z,N);

	surfaceNumber = N - A;
	surfaceIndex = calloc(surfaceNumber,sizeof(int));
	surfaceCood = calloc(surfaceNumber,sizeof(int));

	j = 0;
	for(i=0;i<N;i++)
	{
		if(i < A)
		{
			note[index[i]] = 0;
		}
		else
		{
			surfaceIndex[j] = index[i];
			j ++;
		}
	}

	for(i=0;i<surfaceNumber;i++)
		surfaceCood[i] = i;
	for(i=0;i<surfaceNumber;i++){
		point = rand() % (surfaceNumber-i);
		temp = surfaceCood[i];
		surfaceCood[i] = surfaceCood[point+i];
		surfaceCood[point+i] = temp;
		if(surfaceCood[i] < B)
			surfaceCood[i] = 1;
		else 
			surfaceCood[i] = 2;
	}

	for(i = 0;i < surfaceNumber;i ++)
		note[surfaceIndex[i]] = surfaceCood[i];
	free(index);
}

//相位分离 二合金
void PhaseSeparationNoteInt2(int *note,int N,int A,double *x,double *y,double *z)
{
	int i;
	int *index;

	index = OrderCoodAlongX(x,y,z,N);
	for(i=0;i<N;i++)
	{
		if(i < A)
			note[index[i]] = 0;
		else
			note[index[i]] = 1;
	}
	free(index);
}

//L0排列 二合金
void L0NoteInt2(int *note,int N,int A,double *x,double *y,double *z)
{
	int i;
	int *index;
	double tempX;
	index = OrderCoodAlongX(x,y,z,N);
	for(i=0;i<N;i++)
	{
		tempX = x[index[i]] / 3.924;
		if((tempX-floor(tempX))<=1e-3 && A > 0)
		{
			note[index[i]] = 0;
			A --;
		}
		else
		{
			note[index[i]] = 1;
		}
	}
	free(index);
}

void ShellByShellNoteInt2(int *note, int N, double *x,double *y, double *z)
{
	int i;
	int *shell;
	shell = Shell_Cood(x,y,z,N);
	for(i = 0; i < N; i++)
	{
		note[i] = (shell[i]%2==0)?0:1;
	}
	free(shell);
}

//每个原子的shell
int* Shell_Cood(double *x,double *y,double *z,int N){
	double *R;
	int *neighbour,*flag,*shell;
	int i,j;
	double a0;
	int sh = 0;
	int atom = 0;

	R = calloc(N*N,sizeof(double));
	Distance(x,y,z,R,N);
	a0 = Rmin*sqrt(2);

	neighbour = calloc(N,sizeof(int));
	flag = calloc(N,sizeof(int));
	shell = calloc(N,sizeof(int));

	for(i=0;i<N;i++){
		neighbour[i] = 0;
		flag[i] = 0;
		shell[i] = 0;
	}

	while(1)
	{	
		if(atom == N)
			break;
		for(i=0;i<N;i++)
			neighbour[i]=0;
		for(i=0;i<N-1;i++)
		{	
			for(j=i+1;j<N;j++)
			{	
				if(flag[i]==1||flag[j]==1)
					continue;
				if(*(R+i*N+j)<0.8*a0)
				{
					neighbour[i]++;
					neighbour[j]++;
				}
			}
		}

		for(i=0;i<N;i++)
		{	
			if(flag[i]==1)
				continue;
			if(neighbour[i]<12)
			{
				shell[i] = sh;
				flag[i] = 1;
				atom++;	
			}
		}

		sh++;
	}
	free(R);
	free(flag);
	free(neighbour);
	return shell;
}


int* Shell_Shape(char *shape,int N)
{
	int *shell;
	double *x,*y,*z;

	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));

	ReadCood(shape,N,x,y,z);

	shell = Shell_Cood(x,y,z,N);

	return shell;
} 

//计算距离
void Distance(double *x,double *y,double *z,double *R,int N){
	int i,j;

	Rmax = 0;
	Rmin = 10000;
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);;
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
			if(*(R+i*N+j)<Rmin)
				Rmin = *(R+i*N+j);
			if(*(R+i*N+j)>Rmax)
				Rmax = *(R+i*N+j);
		}
	}
}

void Distance1(COOD *cood,COODDIS* dis)
{
	int i,j;
	int N = cood->N;
	double *R = dis->R;
	double *x = cood->x,*y = cood->y,*z = cood->z;

	dis->N = N;
	dis->Rmax = 0;
	dis->Rmin = 10000;
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			*(R+i*N+j)=(x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j]);;
			*(R+i*N+j)= sqrt(*(R+i*N+j));
			*(R+j*N+i) = *(R+i*N+j);
			if(*(R+i*N+j)<dis->Rmin)
				dis->Rmin = *(R+i*N+j);
			if(*(R+i*N+j)>dis->Rmax)
				dis->Rmax = *(R+i*N+j);
		}
	}
}

char* StoragePath(char *shape,int N,ALLOY *alloy,ATOMNUM *atomNum,char *Output)
{
	int i;
	char Line[20];
	char *Line_End;
	
	//开辟足够大的空间来存储路径
	Line_End = calloc(100,sizeof(char));
	
	//路径: [Output]
	strcpy(Line_End,Output);
	mkdir(Line_End);

	//路径：[Output]\\原子类型
	strcat(Line_End,"\\");
	for( i = 0; i < alloy->atomTypeCount; i++ )
		strcat(Line_End,GetAtomPara(alloy->atoms[i]).name);
	mkdir(Line_End);

	//路径：[Output]\\原子类型\\构型名称
	strcat(Line_End,"\\");
	strcat(Line_End,shape);
	mkdir(Line_End);

	//路径：[Output]\\原子类型\\构型名称\\原子总数
	sprintf(Line,"\\%d",N);
	strcat(Line_End,Line);
	mkdir(Line_End);

	//路径：[Output]\\原子类型\\构型名称\\原子总数\\各原子比例
	strcat(Line_End,"\\");
	for( i = 0; i < atomNum->atomTypeCount; i++ )
	{
		if( i == 0 )
			sprintf(Line,"%.3f",(double)atomNum->numberOfAtom[i] / N);
		else
			sprintf(Line,"-%.3f",(double)atomNum->numberOfAtom[i] / N);
		strcat(Line_End,Line);
	}
	mkdir(Line_End);

	return Line_End;
}

void printDiamond(int *note,int N,COOD *cood,ALLOY *alloy,char *path)
{
	FILE *fp;
	int i;
	
	fp = fopen(path,"w");
	fprintf(fp,"%d\n",N);
	for(i = 0; i < alloy->atomTypeCount; i++)
		fprintf(fp,"%s",GetAtomPara(alloy->atoms[i]).name);
	fprintf(fp,",NCS\n");
	for(i=0;i<N;i++){
		fprintf(fp,"%s\t%lf\t%lf\t%lf\n",GetAtomPara(alloy->atoms[note[i]]).name,cood->x[i],cood->y[i],cood->z[i]);
	}
	fclose(fp);
}


void printResult(int *note,int N,COOD *cood,char *path)
{
	FILE *fp;
	int i;

	fp = fopen(path,"w");
	fprintf(fp,"%d\n",N);
	for(i=0;i<N;i++)
		fprintf(fp,"%d\t%lf\t%lf\t%lf\n",note[i],cood->x[i],cood->y[i],cood->z[i]);	
	fclose(fp);
}



int* orderCoodFromCore(double *x,double *y,double *z,int N)
{
	int i,j;
	double tempR;
	double *R;
	double middleX,middleY,middleZ;
	double maxX,minX,maxY,minY,maxZ,minZ;
	int tempI;
	int *index;
	

	maxX = x[0];maxY = y[0];maxZ = z[0];
	minX = maxX;minY = maxY;minZ = maxZ;

	index = calloc(N,sizeof(int));
	R = calloc(N,sizeof(double));
	for(i=0;i<N;i++)
	{
		maxX = (x[i]>maxX)?x[i]:maxX;
		maxY = (y[i]>maxY)?y[i]:maxY;
		maxZ = (z[i]>maxZ)?z[i]:maxZ;
		minX = (x[i]<minX)?x[i]:minX;
		minY = (y[i]<minY)?y[i]:minY;
		minZ = (z[i]<minZ)?z[i]:minZ;
	}

	middleX = (maxX+minX)/2;
	middleY = (maxY+minY)/2;
	middleZ = (maxZ+minZ)/2;

	for(i = 0;i < N;i ++)
	{
		R[i] = sqrt((x[i]-middleX)*(x[i]-middleX)+(y[i]-middleY)*(y[i]-middleY)+(z[i]-middleZ)*(z[i]-middleZ));
		index[i] = i;
	}

	for(i = 0;i < N - 1;i ++)
		for(j = 0;j < N-1-i;j ++)
		{
			if(R[j] > R[j+1])
			{
				tempR = R[j+1];
				R[j+1] = R[j];
				R[j] = tempR;

				tempR = x[j+1];
				x[j+1] = x[j];
				x[j] = tempR;

				tempR = y[j+1];
				y[j+1] = y[j];
				y[j] = tempR;

				tempR = z[j+1];
				z[j+1] = z[j];
				z[j] = tempR;

				tempI = index[j+1];
				index[j+1] = index[j];
				index[j] = tempI;
			}
		}

	for(i = 0;i < N;i ++)
	{
		R[i] = sqrt((x[i]-middleX)*(x[i]-middleX)+(y[i]-middleY)*(y[i]-middleY)+(z[i]-middleZ)*(z[i]-middleZ));
	}
	free(R);

	return index;
}

//根据x从小到大排列坐标
int*  OrderCoodAlongX(double *x,double *y,double *z,int N)
{
	int i,j;
	int *index;
	double *R;
	double tempR;
	int tempI;
	index = calloc(N,sizeof(int));
	R = calloc(N,sizeof(double));
	for(i = 0;i < N;i ++)
	{
		R[i] = x[i];
		index[i] = i;
	}

	for(i = 0;i < N - 1;i ++)
		for(j = 0;j < N-1-i;j ++)
		{
			if(R[j] > R[j+1])
			{
				tempR = R[j+1];
				R[j+1] = R[j];
				R[j] = tempR;

				tempI = index[j+1];
				index[j+1] = index[j];
				index[j] = tempI;
			}
		}

	free(R);
	return index;
}

//根据R返回各个原子的配位数
int *coordinationNumberR( COODDIS *dis )
{
	int i,j,N;
	int *neighbour;
	double *R;
	double mina0;

	N = dis->N;
	R = dis->R;
	mina0 = dis->Rmin + 0.2;
	neighbour = calloc( N, sizeof(int) );

	for( i = 0; i < N-1; i++ ){	
		for( j = i+1; j < N; j++ ){	
			if( *(R+i*N+j) < mina0 ){
				neighbour[i]++;
				neighbour[j]++;
			}
		}
	}

	return neighbour;
}


void CoordinateNumber( char *input, int N, char *output )
{
	COOD cood;
	COODDIS dis;
	int *note, *neighbour;
	int i,j,sum,typeNumber;
	int num[12] = {0,};
	int (*coord)[12];
	FILE *fp;
	
	Cood_Init( &cood, N );
	note = calloc(N,sizeof(int));
	dis.R = calloc(N*N,sizeof(double));
	
	ReadFile1( input, note, &cood, N );
	
	typeNumber = 0;
	for( i = 0; i < N; i++ )
	{
		if( note[i] > typeNumber )
			typeNumber = note[i];
	}
	typeNumber ++;
	coord = calloc( typeNumber, sizeof(*coord) );
	memset( coord, 0, sizeof(*coord)*typeNumber );

	Distance1( &cood, &dis );
	neighbour = coordinationNumberR( &dis );
	for( i = 0; i < N; i++){
		num[neighbour[i]-1] ++;
		coord[note[i]][neighbour[i]-1] ++;
	}
	
	sum = 0;
	for( i = 0; i < 12; i++ )
		sum = sum + num[i];

	printf( "配位数分析：\n");
	printf( "配位数\tall" );
	for( i = 0; i < typeNumber; i++ )
		printf( "\tatom%d", i );
	printf( "\n" );
	for( i = 0; i < 12; i++ )
	{
		printf( "%d\t%d",i+1,num[i]);
		for( j = 0; j < typeNumber; j++ )
			printf( "\t%d", coord[j][i] );
		printf( "\n" );
	}
	printf("all\t%d\n",sum );

	fp = fopen( output, "w" );
	fprintf( fp, "配位数\tall");
	for( i = 0; i < typeNumber; i++ )
		fprintf( fp, "\tatom%d", i );
	fprintf( fp, "\n");
	for( i = 0; i < 12; i++ )
	{
		fprintf(fp,"%d\t%d", i+1, num[i] );
		for( j = 0; j < typeNumber; j++ )
			fprintf( fp, "\t%d", coord[j][i] );
		fprintf( fp, "\n" );
	}
	fprintf( fp, "all\t%d\n", sum );
	fclose(fp);
	
	free( coord );
	free( neighbour );
	free(note);
	free(dis.R);
	Cood_Free( &cood );
}

void pairCorrelation(char *input,int N,char *Output)
{
	int *note;
	COOD cood;
	double *R;
	double maxX,minX,maxY,minY,maxZ,minZ;
	double middleX,middleY,middleZ;
	int i,midAtom,typeNumber;
	double r,dr;
	int *atom;
	FILE *fp;

	note = calloc(N,sizeof(int));
	Cood_Init( &cood, N );
	R = calloc(N,sizeof(double));
	
	ReadFile1( input, note, &cood, N );
	
	typeNumber = 0;
	for( i = 0; i < N; i++ )
	{
		if( note[i] > typeNumber )
			typeNumber = note[i];
	}
	typeNumber ++;
	atom = calloc( typeNumber, sizeof(int) );
	memset( atom, 0, sizeof(int)*typeNumber );

	maxX = cood.x[0];maxY = cood.y[0];maxZ = cood.z[0];
	minX = maxX;minY = maxY;minZ = maxZ;
	for(i=0;i<N;i++)
	{
		maxX = ( cood.x[i] > maxX )? cood.x[i]:maxX;
		maxY = ( cood.y[i] > maxY )? cood.y[i]:maxY;
		maxZ = ( cood.z[i] > maxZ )? cood.z[i]:maxZ;
		minX = ( cood.x[i] < minX )? cood.x[i]:minX;
		minY = ( cood.y[i] < minY )? cood.y[i]:minY;
		minZ = ( cood.z[i] < minZ )? cood.z[i]:minZ;
	}

	middleX = (maxX+minX)/2;
	middleY = (maxY+minY)/2;
	middleZ = (maxZ+minZ)/2;
	
	midAtom = -1;
	for( i = 0; i < N; i++ )
	{
		if( fabs(cood.x[i]-middleX)<0.001 && fabs(cood.y[i]-middleY)<0.001 && fabs(cood.z[i]-middleZ)<0.001)
		{
			midAtom = i;
			break;
		}
	}
	for(i=0;i<N;i++)
	{
		if( midAtom == -1 )
			R[i] = sqrt((cood.x[i]-middleX)*(cood.x[i]-middleX)+(cood.y[i]-middleY)*(cood.y[i]-middleY)+(cood.z[i]-middleZ)*(cood.z[i]-middleZ));
		else
			R[i] = sqrt((cood.x[i]-cood.x[midAtom])*(cood.x[i]-cood.x[midAtom])+(cood.y[i]-cood.y[midAtom])*(cood.y[i]-cood.y[midAtom])+(cood.z[i]-cood.z[midAtom])*(cood.z[i]-cood.z[midAtom]));
	}
	
	printf("径向分布函数：\n");
	r = 0;
	dr = 0.8;
	fp = fopen( Output, "w" );
	printf( "R" );
	fprintf( fp, "R");
	for( i = 0; i < typeNumber; i++ )
	{
		fprintf( fp, "\tatom%d", i );
		printf( "\tatom%d", i );
	}
	printf( "\n" );
	fprintf( fp, "\n");

	for(r = 0;r < 50;r = r + dr)
	{	
		memset( atom, 0, sizeof(int)*typeNumber );

		for( i = 0;i < N; i++ )
		{
			if(R[i] >= r && R[i] < r + dr)
				atom[note[i]]++;
		}
		printf( "%f", r );
		fprintf( fp, "%f", r );
		for( i = 0; i < typeNumber; i++ )
		{
			printf( "\t%d", atom[i] );
			fprintf( fp, "\t%d", atom[i] );
		}
		printf( "\n" );
		fprintf( fp, "\n" );
	}

	fclose( fp );
	free( atom );
	free( R );
	Cood_Free( &cood );
}


//分层统计
void ShellCount(char *input,int N,char *output)
{
	int i,j,sh,atom,typeNumber,rrr;
	int *neighbour,*flag,*shell,*note;
	int *nums,*r;
	FILE *fp;
	double a0;
	COOD cood;
	COODDIS dis;

	note = calloc(N,sizeof(int));
	Cood_Init( &cood, N );
	dis.R = calloc(N*N,sizeof(double));
	neighbour = calloc(N,sizeof(int));
	flag = calloc(N,sizeof(int));
	shell = calloc(N,sizeof(int));

	ReadFile1( input, note, &cood, N );
	Distance1( &cood, &dis );
	
	typeNumber = 0;
	for( i = 0; i < N; i++ )
	{
		if( note[i] > typeNumber )
			typeNumber = note[i];
	}
	typeNumber ++;
	nums = calloc( typeNumber, sizeof(int) );
	memset( nums, 0, sizeof(int)*typeNumber );
	r = calloc( typeNumber, sizeof(int) );
	memset( r, 0, sizeof(int)*typeNumber );

	a0 = dis.Rmin * sqrt(2);
	for( i = 0; i < N; i++ ){
		neighbour[i] = 0;
		flag[i] = 0;
		shell[i] = 0;
	}
	
	printf( "核壳分布：\n" );
	fp = fopen( output, "w" );
	printf( "shell\tall" );
	fprintf( fp, "shell\tall" );
	for( i = 0; i < typeNumber; i++ )
	{
		printf( "\tatom0%d", i );
		fprintf( fp, "\tatom%d" ,i );		
	}
	printf( "\n" );
	fprintf( fp, "\n" );
	
	atom = 0; 
	sh = 0;
	while(1)
	{	
		if( atom == N ) break;
		
		for( i = 0; i < N; i++ )
			neighbour[i] = 0;
		for( i = 0; i < N-1; i++ )
		{	
			for( j = i+1; j < N; j++ )
			{	
				if( flag[i]==1 || flag[j]==1 )
					continue;
				if( *(dis.R+i*N+j)<0.8*a0 )
				{
					neighbour[i]++;
					neighbour[j]++;
				}
			}
		}

		memset( r, 0, sizeof(int)*typeNumber );
		rrr = 0;
		for( i = 0; i < N; i++ )
		{	
			if( flag[i] == 1 )
				continue;
			if( neighbour[i] < 12 )
			{
				shell[i] = sh;
				flag[i] = 1;
				atom++;	
				r[note[i]] ++;
				rrr ++;
			}
		}

		for( i = 0; i < typeNumber; i++ )
			nums[i] = nums[i] + r[i];
		
		printf( "%d\t%d", sh, rrr );
		fprintf( fp,"%d\t%d", sh, rrr );
		for( i = 0; i < typeNumber; i++ )
		{
			printf( "\t%d", r[i] );
			fprintf( fp,"\t%d", r[i] );		
		}
		printf( "\n" );
		fprintf( fp, "\n" );

		sh++;
	}

	printf( "ALL\t%d", atom );
	fprintf( fp, "ALL\t%d", atom );
	for( i = 0; i < typeNumber; i++ )
	{
		printf( "\t%d", nums[i] );
		fprintf( fp, "\t%d", nums[i] );	
	}
	printf( "\n" );
	fprintf( fp, "\n" );
	fclose(fp);
	
	free( nums );
	free( r );
	free( neighbour );
	free( flag );
	free( shell );
	free( dis.R );
	Cood_Free( &cood );
	free( note );
}

void JAJB(char *input,int N)
{
	int i,j;
	double a0;
	double *R;
	int *note;
	double *x,*y,*z;
	int NAB,NBA,NA,NB;
	int An,Bn;
	double Ja,Jb;
	double Po,Ro,Pr,Rr;
	note = calloc(N,sizeof(int));
	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	R = calloc(N*N,sizeof(double));
	ReadFile(input,note,x,y,z,N);
	Distance(x,y,z,R,N);
	a0 = Rmin * sqrt(2);

	An = 0;Bn = 0;
	for(i=0;i<N;i++)
	{
		if(note[i]==0)
			An++;
		else
			Bn++;
	}
	NAB = 0;NBA = 0;NA = 0;NB = 0;
	for(i=0;i<N-1;i++)
	{
		for(j=i+1;j<N;j++)
		{
			if(*(R+i*N+j) < a0 * 0.8)
			{
				if(note[i] == 0 && note[j] == 0)
				{
					NA = NA + 2;
				}
				if(note[i] == 0 && note[j] == 1)
				{
					NAB = NAB + 1;
					NBA = NBA + 1;
					NA = NA + 1;
					NB = NB + 1;
				}
				if(note[i] == 1 && note[j] == 0)
				{
					NAB = NAB + 1;
					NBA = NBA + 1;
					NA = NA + 1;
					NB = NB + 1;
				}
				if(note[i] == 1 && note[j] == 1)
				{
					NB = NB + 2;
				}
			}
		}// End of for(j=i+1;j<N;j++) 
	}// End of for(i=0;i<N-1;i++) 
	Po = (double)NAB / NA;
	Ro = (double)NBA / NB;
	Pr = (double)An / N;
	Rr = (double)Bn / N;
	Ja = Po/Pr;
	Jb = Ro/Rr;
	printf("Po = %f\tRo = %f\n",Po,Ro);
	printf("Pr = %f\tRr = %f\n",Pr,Rr);
	printf("Ja = %f\tJb = %f\n",Ja,Jb);
}


void printData(char *Line_Date,int N)
{
	char Line_Input[200];
	char Line_Output[200];
	strcpy(Line_Input,Line_Date);
	strcat(Line_Input,"\\result.txt");

	strcpy(Line_Output,Line_Date);
	strcat(Line_Output,"\\CoordinateNumber.xls");
	CoordinateNumber(Line_Input,N,Line_Output);

	strcpy(Line_Output,Line_Date);
	strcat(Line_Output,"\\PairCorrelation.xls");
	pairCorrelation(Line_Input,N,Line_Output);

	strcpy(Line_Output,Line_Date);
	strcat(Line_Output,"\\ShellCount.xls");
	ShellCount(Line_Input,N,Line_Output);

}

/*
double Johnson_F[NN];

void SetJohnsonEnergy()
{
	//FeCu
	double re[2] = {2.481987,2.556162};
	double fe[2] = {1.885957,1.554485};
	double rhoe[2] = {20.041463,22.150141};
	double alpha[2] = {9.81870,7.669911};
	double beta[2] = {5.236411,4.090619};
	double A[2] = {0.392811,0.327584}; 
	double B[2] = {0.64624,0.468735};
	double kappa[2] = {0.170306,0.431307};
	double lambda[2] = {0.340613,0.86214};
	double Fn0[2] = {-2.534992,-2.176490};
	double Fn1[2] = {-0.059605,-0.140035};
	double Fn2[2] = {0.193065,0.285621};
	double Fn3[2] = {-2.282322,-1.750834};
	double F0[2] = {-2.54,-2.19};
	double F1[2] = {0,0};
	double F2[2] = {0.200269,0.702991};
	double F3[2] = {-0.148770,0.683705};
	double eta[2] = {0.391750,0.921150};
	double Fe[2] = {-2.539945,-2.191675};

	int i;
	double XI;


	// 将指数计算做成数组 
	for(i=0;i<NN;i++)
	{
		XI = (double)i/20000;
		NUM1[i] = pow(XI,n1); NUM2[i] = pow(XI,n2); NUM3[i] = pow(XI,n3);
		NUM12[i] = pow(XI,n12); NUM13[i] = pow(XI,n13); NUM23[i] = pow(XI,n23);
		MUM1[i] = pow(XI,m1); MUM2[i] = pow(XI,m2); MUM3[i] = pow(XI,m3);
		MUM12[i] = pow(XI,m12); MUM13[i] = pow(XI,m13); MUM23[i] = pow(XI,m23);
	}
}*/


void saveMatrix(double *x, int N, char* output)
{
	int i,j;
	FILE *f;
	f = fopen(output,"w");
	
	for(i=0;i<N;i++) {
		for(j=0;j<N;j++)
		{	
			if(i>j)
				fprintf(f,"%le\t",*(x + j*N + i));
			else
				fprintf(f,"%le\t",*(x + i*N + j));
		}
		fprintf(f,"\n");
	}

	fclose(f);
}

double getLatticeParameter(ALLOY *alloy)
{
	int i;
	double atomLattice = 0;
	for( i = 0; i < alloy->atomTypeCount; i++)
	{
		atomLattice += GetAtomPara(alloy->atoms[i]).a;
	}
	atomLattice /= alloy->atomTypeCount;
	return atomLattice;
}



//最速下降法弛豫
void SD_File(char *input, ATOM *atoms, int atomTypeCount, int N)
{
	int i=0,j=0;
	int tig = 0;
	int *note;
	double *x,*y,*z,*xx,*yy,*zz;
	double *R,*RR;
	double *FX,*FY,*FZ;
	double E0,E1,F0,F1;
	double eps = 1e-10;
	double nowstep = 0.001;
	int RSTEP = 1000;
	int clocks = 0;
	ALLOY alloy;
	
	alloy.atoms = atoms;
	alloy.atomTypeCount = atomTypeCount;

	note = calloc(N,sizeof(int));
	x = calloc(N,sizeof(double));
	y = calloc(N,sizeof(double));
	z = calloc(N,sizeof(double));
	R = calloc(N*N,sizeof(double));
	FX = calloc(N,sizeof(double));
	FY = calloc(N,sizeof(double));
	FZ = calloc(N,sizeof(double));
	xx = calloc(N,sizeof(double));
	yy = calloc(N,sizeof(double));
	zz = calloc(N,sizeof(double));
	RR = calloc(N*N,sizeof(double));


	ReadFile(input,note,x,y,z,N);
	Distance(x,y,z,R,N);
	E0 = QSCEnergy(note,R,&alloy,N);
	F0 = QSCForce(note,R,x,y,z,FX,FY,FZ,atoms,atomTypeCount,N);
	E1 = E0;
	F1 = F0;
	printf("N = %d\n",N);
	printf("E0 = %lf\n",E0);
	printf("F0 = %lf\n",F0);
	printf("a=%lf\n",Rmin*sqrt(2));


	while(1)
	{   
		//力足够小或迭代次数很大时退出
		if(F0 < eps || clocks >=RSTEP || nowstep < 1e-10){
			break;
		}
		clocks ++;
		//根据力方向移动
		for(i=0;i<N;i++)
		{	
			xx[i] = x[i]+nowstep*FX[i];
			yy[i] = y[i]+nowstep*FY[i];
			zz[i] = z[i]+nowstep*FZ[i];
		}
		
		Distance(xx,yy,zz,RR,N);
		F1 = QSCForce(note,RR,xx,yy,zz,FX,FY,FZ,atoms,atomTypeCount,N);
		E1 = QSCEnergy(note,RR,&alloy,N);
		//如果能量减小接受增大步长，否则步长减小
		if(E1<E0)
		{
			F0 = F1;
			E0 = E1;
			for(i=0;i<N;i++)
			{	
				x[i] = xx[i];
				y[i] = yy[i];
				z[i] = zz[i];
			}
			for(i=0;i<N-1;i++)
			{
				for(j=i+1;j<N;j++)
				{
					*(R+i*N+j) = *(RR+i*N+j);
					*(R+j*N+i) = *(RR+j*N+i);
				}
			}
			nowstep = nowstep *1.05;
			tig = 1;
		}
		else
		{
			nowstep = 0.6*nowstep;
			tig = 0;
		}

       if(tig == 1)
	   {
			printf("clocks=  %d\n",clocks);
			printf("E0    =  %lf\n",E0);
			printf("F     =  %lf\n",F0);
			printf("step  =  %lf\n",nowstep);
		}

	}
	Distance(x,y,z,R,N);
	printf("a=%lf\n",Rmin*sqrt(2));
	//	printDiamond(note,x,y,z,alloy,N,input);

	free(x);
	free(y);
	free(z);
	free(FX);
	free(FY);
	free(FZ);
	free(R);
	free(xx);
	free(yy);
	free(zz);
	free(RR);
}

void Alloy_Init(ALLOY *alloy,...)
{
	va_list argp;
	int argno;
	ATOM temp;

	va_start(argp,alloy);
	argno = 0;
	while( (temp = va_arg(argp,ATOM)) != END )
		argno ++;

	alloy->atomTypeCount = argno;
	alloy->atoms = calloc(argno,sizeof(ATOM));
	va_end(argp);
	
	va_start(argp,alloy);
	argno = 0;
	while( (temp = va_arg(argp,ATOM)) != END )
	{
		alloy->atoms[argno] = temp;
		argno ++;
	}
	va_end(argp);
}

void Alloy_Copy(ALLOY *to, ALLOY *from)
{
	to->atomTypeCount = from->atomTypeCount;
	to->atoms = calloc(to->atomTypeCount,sizeof(ATOM));
	memcpy(to->atoms, from->atoms, to->atomTypeCount * sizeof(ATOM));
}

void Alloy_Free(ALLOY *alloy)
{
	free(alloy->atoms);
	alloy->atoms = NULL;
	alloy->atomTypeCount = 0;
}

void AtomNum_Init(ATOMNUM *atomNum,...)
{
	va_list argp;
	int argno;
	int temp;

	va_start(argp,atomNum);
	argno = 0;
	while( (temp = va_arg(argp,int)) != END )
		argno ++;
	va_end(argp);
	atomNum->atomTypeCount = argno;
	atomNum->numberOfAtom = calloc(argno,sizeof(int));
	
	va_start(argp,atomNum);
	argno = 0;
	while( (temp = va_arg(argp,int)) != END )
	{
		atomNum->numberOfAtom[argno] = temp;
		argno ++;
	}
	va_end(argp);
}

void AtomNum_Copy(ATOMNUM *to, ATOMNUM *from)
{
	to->atomTypeCount = from->atomTypeCount;
	to->numberOfAtom = calloc(to->atomTypeCount,sizeof(int));
	memcpy(to->numberOfAtom, from->numberOfAtom, to->atomTypeCount * sizeof(int));
}

void AtomNum_Free(ATOMNUM *atomNum)
{
	free(atomNum->numberOfAtom);
	atomNum->numberOfAtom = NULL;
	atomNum->atomTypeCount = 0;
}

void Cood_Init(COOD *cood, int N)
{
	cood->x = calloc(N,sizeof(double));
	cood->y = calloc(N,sizeof(double));
	cood->z = calloc(N,sizeof(double));
	cood->N = N;
}

void Cood_Free(COOD *cood)
{
	free(cood->x);
	free(cood->y);
	free(cood->z);
	cood->N = 0;
}









