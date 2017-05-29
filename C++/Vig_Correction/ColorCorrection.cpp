#include<iostream>
#include "math.h"

//corelib
#include "Matrix.h"
#include "ColorCorrection.hpp"
using namespace std;

/* vignetting correction based on radial gradient
*/
int VignettingCorrectionUsingRG(unsigned char* pImage, int ht, int wd, vector<double>& vp)
{
	int halfWd = wd*0.5 ;
	int halfHt = ht*0.5 ;

	int nRadius = sqrt( double(halfHt*halfHt + halfWd*halfWd) ) + 0.5 + 1 ; 

	vp.resize(nRadius);
    
	double* weight = (double*)malloc(ht*wd*sizeof(double));
	//memset(weight, 1, ht*wd*sizeof(double));
	for(int i=0; i<ht*wd; i++)
		weight[i] = 1;

	double* rgImage = (double*)malloc(ht*wd*sizeof(double));
	memset(rgImage, 0, ht*wd*sizeof(double));

	double* A  = (double*)malloc(nRadius*sizeof(double));
	double* At = (double*)malloc(nRadius*sizeof(double));
	double* AtA  = (double*)malloc(nRadius*nRadius*sizeof(double));
	double* sAtA =  (double*)malloc(nRadius*nRadius*sizeof(double));    
	double* sAtL =  (double*)malloc(nRadius*sizeof(double));	
	double* result = (double*)malloc(nRadius*sizeof(double));    
	double* mB = (double*)malloc(nRadius*nRadius*sizeof(double));
	memset(mB, 0, nRadius*nRadius*sizeof(double));
	double* mBt = (double*)malloc(nRadius*nRadius*sizeof(double));
	memset(mBt, 0, nRadius*nRadius*sizeof(double));
	double* mBtB = (double*)malloc(nRadius*nRadius*sizeof(double));
	memset(mBtB, 0, nRadius*nRadius*sizeof(double));
   
	//smooth constrait
	//lambda*(2*numPixels/numR)
	double lamdaS = 0.15*2*(wd*ht)/double(nRadius);
	for(int i=1; i<nRadius-1; i++)
	{
		mB[i*nRadius+i]   = -2;//*lamdaS;
		mB[i*nRadius+i-1] = 1; //*lamdaS;
		mB[i*nRadius+i+1] = 1; //*lamdaS;
	}
	//matrix_print(nRadius, nRadius, mB);
	//PrintMatrix(mB, nRadius, nRadius);
    transpose(mB, mBt, nRadius, nRadius);
	mult(mBt, mB, mBtB, nRadius, nRadius, nRadius);
	//PrintMatrix(mBtB, nRadius, nRadius);
    

	//calculate the radial gradients of image
	double shift = 1;
	double eps = 0.000001;
	for(int j=1; j<ht; j++)
	{
		for(int i=1; i<wd; i++)
		{
			int cx = i - halfWd;
			int cy = j - halfHt;

			//calculate the radius
			int radius = sqrt( double(cx*cx + cy*cy ) ) + 0.5;
			
			//calculate the gradient
			double dx = log( double(pImage[j*wd+i]) + shift ) - log( double(pImage[j*wd+i-1]) + shift );
			double dy = log( double(pImage[j*wd+i]) + shift ) - log( double(pImage[(j-1)*wd+i]) + shift );
			
			//calculate the radial gradient
			double rg = double(cx*dx+cy*dy) / sqrt( double(cx*cx+cy*cy) + eps );
			rgImage[j*wd+i] = rg;
		}
	}
	//PrintMatrix(rgImage, ht, wd);

	//weighted least square solution
	printf("\n WLS... \n");
	for(int iterIndex=0; iterIndex<5; iterIndex++)
	{
		memset(sAtA, 0, nRadius*nRadius*sizeof(double));
		memset(sAtL, 0, nRadius*sizeof(double));

		for(int j=1; j<ht; j++)
			for(int i=1; i<wd; i++)
			{
				//printf("%d %d \n", j, i);

				memset(A,   0, nRadius*sizeof(double));
				memset(At,  0, nRadius*sizeof(double));
				memset(AtA, 0, nRadius*nRadius*sizeof(double));

				//calculate the radius
				int cx = i - halfWd;
				int cy = j - halfHt;
				int radius = sqrt( double(cx*cx + cy*cy ) ) + 0.5;

				double rg = rgImage[j*wd+i];

				//calculate the AtA of each pixel
				double right=0;
				if(radius>0 && radius<nRadius)
				{
					A[radius]   = 1;
					A[radius-1] = -1;
					right = rg;
				}	
				//for (int i = 0; i < nRadius; i++)
				//	cout << A[i] << " ";
				//cout << endl;
				for(int k=0; k<nRadius; k++)
					At[k] = A[k];

				mult(At, A, AtA, nRadius, 1, nRadius);
				
				//PrintMatrix(AtA, nRadius, nRadius);

				//sum of AtA
				double w2 = weight[j*wd+i]*weight[j*wd+i];
				//printf("%d %d %lf \n", j, i, w2);
				for(int k=0; k<nRadius*nRadius; k++)
				{
					sAtA[k] += AtA[k]*w2;
				}
				for(int k=0; k<nRadius; k++)
				{
					sAtL[k] += At[k]*right*w2;
				}
			}

		/////////////////////////////  adding constraints ///////////////////////
		//smooth constraint
		for(int i=0; i<nRadius*nRadius; i++)
			sAtA[i] += lamdaS*lamdaS*mBtB[i];

		//vignetting value constraint, make them close to 1
		double eps = 0.03;
		for(int i=0; i<nRadius; i++)
			sAtA[i*nRadius+i] += eps;
		//////////////////////////////////////////////////////////////////////////
		
		invers_matrix(sAtA, nRadius);
		mult(sAtA, sAtL, result, nRadius, nRadius, 1);

		
		//printf("vignetting paras..... \n");
		for(int i=0; i<nRadius; i++)
		{
			vp[i] = result[i];
			//cout << result[i] << endl;
			//printf("%lf ", vp[i]);
		}
		//update weight
		double alpha = 0.6;
		for(int j=1; j<ht; j++)
			for(int i=1; i<wd; i++)
			{
				double rz; //radial gradient of image
				double rv; //radial gradient of vignetting paras

				int cx = i - halfWd;
				int cy = j - halfHt;
				//calculate the radius
				int radius = sqrt( double(cx*cx + cy*cy ) ) + 0.5;
				radius=max( 1, min(nRadius-1, radius) );

				//rv = log(vp[radius])-log(vp[radius-1]);
				rv = vp[radius] - vp[radius-1];
				rz = rgImage[j*wd+i];

				double s1 = fabs(rz-rv); //sqrt( (rz-rv)*(rz-rv) );
				double s2 = alpha*pow(s1, alpha-1);
				weight[j*wd+i] = exp(-s1)*(1-exp(-s2));
			}
	}    

	//for (int i = 0; i < ht*wd; i++)
	//	cout << weight[i] << endl;
	free(A);
	free(At);
	free(AtA);
	free(sAtA);
	free(sAtL);
	free(mBt);
	free(mB);
	free(mBtB);
	free(result);
	free(weight);
	free(rgImage);

	return 0;
}




