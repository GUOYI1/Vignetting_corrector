
#ifndef CORELIB_COMMONFILE_FUNCTION
#define CORELIB_COMMONFILE_FUNCTION

#include "commondata.h"

//Common functions
# define        sign(A,B)        ((A)<0?-(B):(B))
# define        sgn(A)        ((A)<0?-1:1)

#define CHARN 256


//2D arrays
double **f2d (int nr, int nc);
int **f2i (int nr, int nc);
BYTE **f2b (int nr, int nc);
char **f2c (int nr, int nc);
MYCOLOR **f2m(int nr, int nc); 
short **f2s (int nr, int nc);
float **f2f(int nr, int nc);


//Release memory
void FreeArray_double(double **array, int nr, int nc);
void FreeArray_int(int **array, int nr, int nc);
void FreeArray_short(short **array, int nr, int nc);
void FreeArray_BYTE(BYTE **array, int nr, int nc);
void FreeArray_COLOR(MYCOLOR **array, int nr, int nc);
void FreeArray_char(char **array, int nr, int nc);
void FreeArray_float(float **array, int nr, int nc);

//Smooth the pixels of image
double Horn_t(BYTE **image, int y, int x, int ht, int wd);
int iHorn_t(BYTE **image, int y, int x, int ht, int wd);
double Horn_t(int **image, int y, int x, int ht, int wd);
int iHorn_t(int **image, int y, int x, int ht, int wd);
double Horn_t3(MYCOLOR **image, int y, int x, int ht, int wd, int pixel);
int iHorn_t3(MYCOLOR **image, int y, int x, int ht, int wd, int pixel);

//Calculate the Nth power of a number
double CalPower(double x, int n);

//Delete the files under a certain path
void DeleteFileInPath(char* pathname);
void DeleteAllFileInPath(char* pathname);

//Judge whether a point is in a rectangle
bool PtInRect(const MyRect *lprc, int x, int y);

//Get the distance between 2 pixels
int GetColorPixelDis(MYCOLOR pix1, MYCOLOR pix2);
//Get the Next point on the spiral
MyPoint GetNextSpiralPoint(MyPoint inpt, int dis);
//Round function
//int round(double x);

void GetImageFileName(char *filename, char *inpath, int &n, int &nfile, char type, bool flag);
void GetDirFileName( char** filename, char *inpath, int *n, int *nfile, const char* type, bool flag);

//get all files in directory, including all subdirectories
void GetDirFileNameAll( char** filename, char *inpath, int *n);


int  MyGetFileSize(char* filename);


//for sfm project progress bar 
int ReadProgressValueFile(double& dValue);
int WriteProgressValueToFile(double nValue);

/*
     d:\data\1.jpg ---> d:\data\ + "subpath\" + 1. + "postfix"
*/
int GenerateProductFile(char* srcfile, char* subpath, char* postfix, char** productFile);



/*
     d:\data\1.jpg ---> d:\data\ + "subpath"
*/
int GenerateProductPath(char* srcfile, char* subpath, char** productPath);


/*
     d:\data\1.jpg ---> 1
*/
int GetTitleName(char* srcfile, char** titleName);




/*
     d:\data\1.jpg ---> 1.dat
*/
int GetTitleWithPosfix(char* srcfile, char** titleName);



#endif