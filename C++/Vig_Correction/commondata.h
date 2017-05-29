

#ifndef COMMONDATA_H
#define COMMONDATA_H 

#include <vector>
using namespace std;


#define MIN_DATA 1e-10
#define MAX_CURVE_POINT 1024

typedef unsigned char byte;
typedef unsigned char BYTE;


#ifndef NULL
#define NULL 0
#endif

//图像尺寸参数
//#define U9_IMG_WIDTH   320  //处理图像的宽度，固定为320，部分代码需要乘以320，可改为两次移位实现。
//#define U8_IMG_HEIGHT  240  //处理图像的高度，固定位240
#define U17_IMG_SIZE   76800  //处理图像的象素点数，等于宽度和高度的乘积

//前景矩形搜索参数
#define U6_MAX_ROW_EDGE_NUM   1024   //二值图像每行包含的最大连续段数
#define U13_MAX_EDGE_NUM      76800   //一幅图像中包含的最大连续段数（＝IMG_HEIGH*MAX_ROW_EDGE_NUM）
#define U8_MAX_LABEL_SUM      100000  //二值面具中前景矩形的最大数量 

#define  U8_MAXMEASURE        100000  //返回到物体跟踪模块的前景矩形的最大数量

//#define U13_MAX_EDGE_NUM      10000


typedef struct stMySize
{
	int height, width;
}MySize;


typedef struct stMyPointF
{
	double x,y;	
}MyPointF;

typedef struct 
{
	int x;
	int y;
}iPOINT;

typedef struct stMyPoint
{
	int x,y;
}stPoint,iPOINT2,MyPoint;

typedef struct stMyPoint3d
{
	unsigned short x,y,z;
}MyPoint3d;

typedef struct stPOINT3
{
	double x,y,z;
	bool   bIsOutSide;
	double f;
}POINT3, POINT3D;


/*
typedef struct stPOINT3d
{
	double x,y,z;
}POINT3D;
*/

/*
typedef struct stPoint2D
{
	double x,y;
}Point2D;
*/


typedef struct stPOINT2	
{
	double x;
	double y;
	int   type;
	int   code;
	int   mask;
}POINT2;

/*
typedef struct stEDGE1D
{
	int left;
	int right;
	int type;
}EDGE1D;
*/

//前景段参数EDGE1D
//为3.3节前景矩形搜索模块的前景段定义的结果，前两个为前景段左侧和右侧点的x坐标，u8_index为该段的标识。
//这个标识不能超过前面定义的U8_MAX_LABEL_SUM＝128。这个结构总的位数为26。
typedef struct EDGE1D
{
	unsigned short left;
	unsigned short right;
	unsigned short u8_index;//可改为6位实现，同时将U8_MAX_LABEL_SUM改为64
	int index;
}EDGE1D; //请注意这个结构可以缩小为24位，即将u8_index改为u6_index，U8_MAX_LABEL_SUM改为64，
		 //这样对大部分视频均可工作，不过当视频非常复杂时会出些问题


typedef struct stCOLOR
{
   unsigned char r,g,b;
}COLOR,MYCOLOR;

typedef struct stRECT
{
	//int l,t,r,b;
	unsigned short left, top, right, bottom;
	int x,y;
	int ht,wd;
	int height, width;
}MyRect;

typedef struct stRECTF
{
	double x, y, width, height;
	double left, top, right, bottom;
}MyRectF;


typedef struct stLabelRect 
{
	unsigned short left;
	unsigned short top;
	unsigned short right;
	unsigned short bottom;
	unsigned char  bIsValid;
}LabelRect;

//stucture for Marching cube
typedef struct stMCPt
{
	float x,y,z;
	bool  bIsOutSide;
	float f;
}MCPt;


template<class T>
struct PYRAMID_IMAGE_GENERAL
{
	int   nL;
	T** pLevelImage;       //pointer to pyramid images of all leves
	unsigned char** pMask; //for mosaic
	int*  pHt;
	int*  pWd;
};


struct IMAGE_GENERAL
{
	short*   pImage;     //image pointer
	int      nHt;
	int      nWd;
};



//typedef unsigned char BOOL;
typedef unsigned char BYTE;
#define false 0
#define true  1

  
#define max(a,b)    (((a) > (b)) ? (a) : (b))
#define min(a,b)    (((a) < (b)) ? (a) : (b))
#define WIDTHBYTES(bits)        ((((bits) + 31)>>5)<<2)

#define MAX_VALUE   100000000
#define PI  3.1415926
#define DPI 0.01745329
//#define DPI 0.017453292519943295769236907684886

#define MAX_COUNT 2000

//affine similarity macro
#define SIMILARITY_X(x,y,p) p[0]*x-p[1]*y+p[2]
#define SIMILARITY_Y(x,y,p) p[1]*x+p[0]*y+p[3]

//affine transform macro
#define AFFINE_X(x,y,p) p[0]*x+p[1]*y+p[2]
#define AFFINE_Y(x,y,p) p[3]*x+p[4]*y+p[5]

#define POLYNOMIAL_X(x,y,p) p[0]*x + p[1]*y + p[2]*x*y + p[3]*x*x + p[4]*y*y + p[5]
#define POLYNOMIAL_Y(x,y,p) p[6]*x + p[7]*y + p[8]*x*y + p[9]*x*x + p[10]*y*y + p[11]

#define HOMOGRAPHY_X(x,y,H) (H[0]*x+H[1]*y+H[2]) / ( H[6]*x+H[7]*y+H[8] )
#define HOMOGRAPHY_Y(x,y,H) (H[3]*x+H[4]*y+H[5]) / ( H[6]*x+H[7]*y+H[8] )

#define AngleToRadian(angle)   angle/180.0*PI
#define RadianToAngle(radian)  radian/PI*180.0

//#define RECT_WD(rect)    

/* Some math libraries may not have the sincos() routine */
#ifdef _SINCOS_
//void sincos();
#define SINCOS(x,s,c)   sincos(x,&s,&c)
#else
//double sin(), cos();
#define SINCOS(x,s,c)   s=sin(x);c=cos(x)
#endif


typedef struct structLINE
{
	MyPointF v1,v2;  //vertex of line
	double   rou,sita; //line polar parameters
	double   a,b,c;    //ax+by+c = 0;
	double   len;
}stLINE;


typedef struct structEllipse
{
	double x0, y0; //center point
	double a, b;   //the axis length of long and short
	double angle;  //
}stEllipse;


//for image 
typedef  vector<int> ImageMatInt;
typedef  vector<unsigned char> ImageMatByte;



#endif