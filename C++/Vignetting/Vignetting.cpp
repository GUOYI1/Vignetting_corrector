#include<cv.h>
#include<cvaux.h>
#include<cxcore.h>
#include<highgui.h>
#include<math.h>
#include<iostream>

using namespace std;

int main(int argc, char** argv)
{
	IplImage* img = cvLoadImage("RS4.png", -1);
	IplImage* vignette = cvCreateImage(cvSize(img->width, img->height), img->depth, img->nChannels);
	IplImage* V = cvCreateImage(cvSize(img->width, img->height), img->depth, 1);
	uchar* V_Initial = new uchar[img->widthStep*img->height];
	vignette = cvCloneImage(img);
	cvSetData(V,V_Initial,img->widthStep);
	double cx = 0.5*img->width;
	double cy = 0.5*img->height;
	double maxDist = 1 / sqrt(cx*cx+cy*cy);
	uchar* p = NULL;
	uchar* q = NULL;
	for (int i = 0; i < vignette->height; i++)
	{
		p = (uchar*)(vignette->imageData + i*vignette->widthStep);
		q = (uchar*)(V->imageData + i*V->widthStep);
		for (int j = 0; j < vignette->width; j++)
		{
			double dist = sqrt((i - cy)*(i - cy) + (j - cx)*(j - cx));
			double ratio= 0.6 / (1.0 + exp((dist * maxDist - 0.8) * 10)) + 0.4;;
			int n = vignette->nChannels;
			if (n == 1)
				p[j] = p[j] * ratio;
			else
			{
				p[n*j] = p[n*j] * ratio;
				p[n*j+1] = p[n*j+1] * ratio;
				p[n*j+2] = p[n*j+2] * ratio;
			}
			q[j] = 255*ratio;
		}

	}
	cvNamedWindow("Image", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Vignette", CV_WINDOW_AUTOSIZE);
	cvNamedWindow("Estimate", CV_WINDOW_AUTOSIZE);
	cvShowImage("Vignette", vignette);
	cvShowImage("Image", img);
	cvShowImage("Estimate", V);

	cvWaitKey(0);
	cvSaveImage("white4_vig.png", vignette);
	cvSaveImage("white4_Est_truth.png", V);
	cvReleaseImage(&img);
	cvReleaseImage(&vignette);
	delete(V_Initial);
//	cvReleaseData(V);
//	cvReleaseImage(&V);
	cvDestroyAllWindows;
}