//% function Z=poly2d(sz,D,a)
//% 
//% Z=poly2d(sz,D,a) does 2-variat polynomial interpolation given the
//% parameters of the polynmomial, with equation (14) in the journal version
//% of bias correction
//% 
//% Input arguments:
//% 
//% sz: 2-component vector specifying the size of the image, e.g. [254 360]
//% D: scalar value shwoing the degree of the polynormial
//% a: vector containing the polynomial parameters in the order determined by
//%    equation (14) in the journal version of NICESIGN  
//% 
//% Ouput argument:
//% 
//% Z: the interpolation image.
//% 
//% Note:
//% 1. x,y value of the polynormial are 1,2,....,sz(2), 1,2,...,sz(1), repectively
//% 2. this function will call the mex file in name poly2d.dll (windows) etc.
//% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
//% Designed and implemented by Yuanjie Zheng at picsl of UPenn
//% Sep. 07, 2009. contact: zheng.vision@gmail.com
//% All rights reserved
//% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

#include "mex.h"
#include <cmath>


void getSize(const mxArray * pvector,int & sz1, int & sz2);
int getD(const mxArray * pscal);
int geta(const mxArray * pavec,double * & pa);
//Z=poly2d(sz,D,a);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	//parse input parameters
	if(nrhs!=3){
		mexErrMsgTxt("There are totally 3 input paramters for getVector: sz,D,a");
	}
	if (nlhs!=1){
		mexErrMsgTxt("There is only 1 output parameter which is the image in the specified size!");
	}
	//get image size
	int sz1,sz2;
	getSize(prhs[0],sz1,sz2);

	//create output matrix
	plhs[0] = mxCreateDoubleMatrix(sz1,sz2,mxREAL);
	double *pData=mxGetPr(plhs[0]);

	//get D
	int D=getD(prhs[1]);

	//get a
	double *pa=0; int num;
	num=geta(prhs[2],pa);
	if (num!=(D+1)*(D+2)/2)
		mexErrMsgTxt("The length of the input a vecgtor is not correct!");

    //interpolation
	int t,l; double sum; int ind;
	for(int i=0;i<sz1;i++)
		for(int j=0;j<sz2;j++)
		{
			sum=0; ind=0;
			for(t=0;t<=D;t++)
				for(l=0;l<=(D-t);l++){
					sum+=pa[ind]*pow((double)(j+1),t)*pow((double)(i+1),l);
					ind++;
				}
			
			pData[j*sz1+i]=sum;
		}
	//release variables
	delete pa;
}

void getSize(const mxArray * pvector, int & sz1, int & sz2){
		int nDims=mxGetNumberOfDimensions(pvector);
		if(nDims>2)
		{
			mexErrMsgTxt("The input sz for poly2d must be a vector!");
		}
		const int* dims=mxGetDimensions(pvector);;

		if(dims[0]!=1 && dims[1]!=1)
			mexErrMsgTxt("The input sz of poly2d must be a vector!");
		int nPara=dims[0]+dims[1]-1;

        double* para=(double *)mxGetPr(pvector);
		sz1=(int)para[0]; sz2=(int)para[1];
}
int getD(const mxArray * pscal)
{
		int nDims=mxGetNumberOfDimensions(pscal);
		if(nDims>2)
		{
			mexErrMsgTxt("The input D for poly2d must be a scalar!");
		}
		const int* dims=mxGetDimensions(pscal);

		if(!(dims[0]==1 && dims[1]==1))
			mexErrMsgTxt("The input D of poly2d must be a scalar!");

		double* para=(double *)mxGetPr(pscal);
		return (int)para[0];
}
int geta(const mxArray * pavec,double * & pa)
{
		int nDims=mxGetNumberOfDimensions(pavec);
		if(nDims>2)
		{
			mexErrMsgTxt("The input a for poly2d must be a vector!");
		}
		const int* dims=mxGetDimensions(pavec);;

		if(dims[0]!=1 && dims[1]!=1)
			mexErrMsgTxt("The input a of poly2d must be a vector!");
		//length of the a vector
		int nPara=dims[0]+dims[1]-1;

		//apply for space for a
		pa=new double[nPara];
		
		//assign value
		double* para=(double *)mxGetPr(pavec);
		for (int i=0;i<nPara;i++)
			pa[i]=para[i];
		return nPara;
}