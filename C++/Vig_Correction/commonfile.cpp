
//#include "stdafx.h"
#include <Math.h>
#include <io.h>
#include <stdio.h>
#include "string.h"
#include "malloc.h"

#include "CommonData.h"
#include "commonfile.h"


double **f2d (int nr, int nc)
{
	double **x;
	int i;

	x = (double **)calloc ( nr, sizeof (double *) );
	if (x == NULL)	
		return NULL;
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (double *) calloc ( nc, sizeof (double)  );
		if (x[i] == NULL)
			  return NULL;				  
	}
	return x;
}

void FreeArray_double(double **array, int row, int col)
{
	int i;

	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}

int **f2i (int nr, int nc)
{
	int **x;
	int i;

	x = (int **)calloc ( nr, sizeof (int *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (int *) calloc ( nc, sizeof (int)  );
		if (x[i] == NULL)
			  return NULL;				  
	}
	return x;
}

short **f2s (int nr, int nc)
{
	short **x;
	int i;

	x = (short **)calloc ( nr, sizeof (short *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (short *) calloc ( nc, sizeof (short)  );
		if (x[i] == NULL)
			  return NULL;				  
	}
	return x;
}

void FreeArray_int(int **array, int row, int col)
{
	int i;

	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}

void FreeArray_short(short **array, int row, int col)
{
	int i;

	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}

BYTE **f2b (int nr, int nc)
{
	BYTE **x;
	int i;

	x = (BYTE **)calloc ( nr, sizeof (BYTE *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (BYTE *) calloc ( nc, sizeof (BYTE)  );
		if (x[i] == NULL)
			  return NULL;				  
	}

	return x;
}

MYCOLOR **f2m (int nr, int nc)
{
	MYCOLOR **x;
	int i;

	x = (MYCOLOR **)calloc ( nr, sizeof (MYCOLOR *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (MYCOLOR *) calloc ( nc, sizeof (MYCOLOR)  );
		if (x[i] == NULL)
			  return NULL;				  
	}

	return x;
}

float **f2f(int nr, int nc)
{
	float **x;
	int i;

	x = (float **)calloc ( nr, sizeof (float *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (float *) calloc ( nc, sizeof (float)  );
		if (x[i] == NULL)
			  return NULL;				  
	}
	return x;
}

void FreeArray_COLOR(MYCOLOR **array, int row, int col)
{
	int i;

	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}


void FreeArray_BYTE(BYTE **array, int row, int col)
{
	int i;

	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}

char **f2c (int nr, int nc)
{
	char **x;
	int i;

	x = (char **)calloc ( nr, sizeof (BYTE *) );
	if (x == NULL)	
		return NULL;		
	
	for (i=0; i<nr; i++)
	{  
		x[i] = (char *) calloc ( nc, sizeof (BYTE)  );
		if (x[i] == NULL)
			  return NULL;				  
	}
	return x;
}

void FreeArray_char(char **array, int row, int col)
{
	int i;
	for(i = 0;i < row; ++i)
		free(array[i]);
	free(array);
}
void FreeArray_float(float **array, int nr, int nc)
{	
	int i;
	for(i = 0;i < nr; ++i)
		free(array[i]);
	free(array);
}

double CalPower(double x, int n)
{
	int i;
	double res = 1.0;
	for(i = 0;i < n; ++i)
		res *= x;
	return res;
}
/*int round(double x)
{
	return ((int)((x)>0 ? (x)+0.5 : (x)-0.5));	
}

int  MyGetFileSize(char* filename)
{
	int filelen;  
	FILE *fp;
	fp = fopen(filename, "rb");
	fseek(fp, 0L, SEEK_END);  
	filelen = ftell(fp);
	rewind(fp);
	fclose(fp);

	return filelen;
}



void GetDirFileNameAll( char** filename, char *inpath, int *n)
{
	char m_strPathName[256];
	char newfile[256], *filetype;
	_finddata_t findinfo;
	long nFileHandle,nTemp;
	
	strcpy(m_strPathName, inpath);
	strcat(m_strPathName, "\\*.*");	
	//strcat(m_strPathName, type);
	
	nFileHandle = _findfirst(m_strPathName, &findinfo);
	if(nFileHandle == -1)
	{
		printf("Bitmap sequence open error!\n");
		return;
	}
	
	do 
	{					
		if(strcmp(findinfo.name,".")&&strcmp(findinfo.name,".."))
		{
			strcpy(newfile, inpath);
			strcat(newfile, "\\");
			strcat(newfile, findinfo.name);
			if(findinfo.attrib&0x10) //judge the directory
				GetDirFileNameAll(filename, newfile, n);
			else
			{
				strcpy( filename[*n], newfile);
				(*n) ++;							
			}
		}		
	}while(!(nTemp = _findnext(nFileHandle,&findinfo)));		
	_findclose(nFileHandle);
}



int WriteProgressValueToFile(double dValue)
{
	static double dValueSum = 0;

	dValueSum += dValue;
	
	printf("progress value: %lf  %lf\n", dValueSum, dValue);

	char* logfile="c:\\sfm\\progress.txt";
	FILE* fp = fopen(logfile, "w");
	if(fp!=NULL)
		fprintf(fp, "%lf", dValueSum);
	fclose(fp);
	return 0;
}

int ReadProgressValueFile(double& dValue)
{
	char* logfile="c:\\sfm\\progress.txt";
	FILE* fp = fopen(logfile, "r");
	if(fp!=NULL)
		fscanf(fp, "%lf", &dValue);
	fclose(fp);
	return 0;
}


int GenerateProductPath(char* srcfile, char* subpath, char** productPath)
{

	//char productFile[256];
	*productPath = (char*)malloc(512);
	memset(*productPath, '\0', 512);

	char* pdes = strrchr(srcfile, '\\');
	int   index = pdes - srcfile;

	//get the imagepath
	char imagepath[512];
	memset(imagepath, '\0', 512);
	strncpy(imagepath, srcfile, index);

	//get the file title
	char title[512];
	strcpy(title, srcfile+index+1);

	//generate new file
	sprintf(*productPath, "%s\\%s", imagepath, subpath);

	return 0;
}


int GenerateProductFile(char* srcfile, char* subpath, char* postfix, char** productFile)
{
	//char productFile[256];
	*productFile = (char*)malloc(512);
	memset(*productFile, '\0', 512);

	char* pdes = strrchr(srcfile, '\\');
	int   index = pdes - srcfile;
	
	//get the imagepath
	char imagepath[512];
	memset(imagepath, '\0', 512);
	strncpy(imagepath, srcfile, index);

	//get the file title
	char title[512];
	strcpy(title, srcfile+index+1);
	
	//generate new file
	sprintf(*productFile, "%s\\%s\\%s", imagepath, subpath, title);
	strcpy(*productFile+strlen(*productFile)-3, postfix);	

	return 0;
}


int GetTitleName(char* srcfile, char** titleName)
{
	//char productFile[256];
	*titleName = (char*)malloc(512);
	memset(*titleName, '\0', 512);

	char* pdes = strrchr(srcfile, '\\');
	int   index = pdes - srcfile;

	//get the imagepath
	//char imagepath[512];
	//memset(imagepath, '\0', 512);
	//strncpy(imagepath, srcfile, index);

	//get the file title
	//char title[512];
	strcpy(*titleName, srcfile+index+1);
    
	int nLen = strlen(*titleName)-3;

	strcpy( (*titleName)+nLen, '\0');	


	return 0;
}


int GetTitleWithPosfix(char* srcfile, char** titleName)
{
	*titleName = (char*)malloc(512);
	memset(*titleName, '\0', 512);

	char* pdes = strrchr(srcfile, '\\');
	int   index = pdes - srcfile;

	strcpy(*titleName, srcfile+index+1);

	return 0;
}*/