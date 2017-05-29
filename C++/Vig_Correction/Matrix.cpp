
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "matrix.h"

#include "commonfile.h"


/*Matrix  Multiplication*/
void mult(double *m1,double *m2,double *result,int i_1,int j_12,int j_2)
{
	int i,j,k;
	for(i=0;i<i_1;i++)
		for(j=0;j<j_2;j++)
		{
			result[i*j_2+j]=0.0;

			for(k=0;k<j_12;k++)
				result[i*j_2+j] += m1[i*j_12+k]*m2[j+k*j_2];
		}
	return;
}

/*Matrix Inverse*/
int invers_matrix(double *m1,int n)
{ 
	int *is,*js;
	int i,j,k,l,u,v;
	double temp,max_v;
	is=(int *)malloc(n*sizeof(int));
	js=(int *)malloc(n*sizeof(int));
	if(is==NULL||js==NULL)
	{
		printf("out of memory!\n");
		return(0);
	}
	for(k=0;k<n;k++)
	{
		max_v=0.0;
		for(i=k;i<n;i++)
			for (j = k; j < n; j++)
			{
				temp = fabs(m1[i*n + j]);
				if (temp > max_v)
				{
					max_v = temp; is[k] = i; js[k] = j;
				}
			}
		if(max_v==0.0)
		{
			free(is); free(js);
			printf("invers is not availble!\n");
			return(0);
		}
		if(is[k]!=k)
		{
			for (j = 0; j<n; j++) 
			{
				u = k*n + j; v = is[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		}

		if(js[k]!=k)
		{
			for (i = 0; i<n; i++) {
				u = i*n + k; v = i*n + js[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		}

		l=k*n+k;
		m1[l]=1.0/m1[l];
		for(j=0;j<n;j++)
		if(j!=k)
		{
			u=k*n+j;
			m1[u]*=m1[l];
		}
		for (i = 0; i < n; i++)
		{
			if (i != k)
				for (j = 0; j < n; j++)
					if (j != k)
					{
						u = i*n + j;
						m1[u] -= m1[i*n + k] * m1[k*n + j];
					}
		}
		for(i=0;i<n;i++)
		{
			if (i != k)
			{
				u = i*n + k;
				m1[u] *= -m1[l];
			}
		}

	}
	for(k=n-1;k>=0;k--)
	{
		if(js[k]!=k)
		{
			for (j = 0; j<n; j++) {
				u = k*n + j; v = js[k] * n + j;
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		}

		if(is[k]!=k)
		{
			for (i = 0; i<n; i++) {
				u = i*n + k; v = i*n + is[k];
				temp = m1[u]; m1[u] = m1[v]; m1[v] = temp;
			}
		}

	}
	free(is); 
	free(js);
	return(1);
}
/*Matrix Transpose*/
void transpose(double *m1,double *m2,int m,int n)
{ 
	int i,j;
	for(i=0;i<m;i++)
		for(j=0;j<n;j++)
			m2[j*m+i]=m1[i*n+j];
	return;
}


/* Function: calculate Eigenvectors and Eigenvalues sorted by decreasing eigenvaules
*/
void JacobiFeatValue(float **inMatrix, float *featl, float **featv, int N)
{
	int i,j;
	float **S, **R; /*Eigenvector matrix*/
	float maxa,t;
	float mine = MIN_DATA;
	int p, q, sw;
	float omega;
	float m, n;
	float sint2, sint, cost, cost2;
	int *swap;


	swap = (int *)malloc(sizeof(int)*N);

	S = f2f(N, N);
	R = f2f(N, N);
	for(i = 0;i < N; ++i)
	{
		for(j = 0;j < N; ++j)
		{
			S[i][j] = 0;
		}
		S[i][i] = 1;
	}
	while(1)
	{		
		maxa = fabs(inMatrix[0][1]);
		p = 0;q = 1;
		for(i = 0;i < N; ++i)
		{
			for(j = i + 1;j < N; ++j)
			{
				if(fabs(inMatrix[i][j]) > maxa)
				{
					maxa = fabs(inMatrix[i][j]);
					p = i;q = j;
				}
			}
		}

		if(maxa < mine)
			break;
		m = -inMatrix[p][q];
		n = (inMatrix[q][q] - inMatrix[p][p])/2;
		omega = sgn(n)*m/(sqrt(m*m + n*n));
		sint2 = omega;
		cost2 = sqrt(1 - sint2*sint2);
		sint = omega/(sqrt(2*(1 + cost2)));
		cost = sqrt(1 - sint*sint);
		for(i = 0;i < N; ++i)
			memcpy(R[i], inMatrix[i], sizeof(float)*N);
		R[p][p] = inMatrix[p][p]*cost*cost + inMatrix[q][q]*sint*sint + inMatrix[p][q]*sint2;
		R[q][q] = inMatrix[p][p]*sint*sint + inMatrix[q][q]*cost*cost - inMatrix[p][q]*sint2;
		R[q][p] = R[p][q] = (inMatrix[q][q] - inMatrix[p][p])*sint2/2 + inMatrix[p][q]*cost2;
		for(i = 0;i < N; ++i)
		{
			if(i == q||i == p)
				continue;
			R[i][p] = R[p][i] = inMatrix[p][i]*cost + inMatrix[q][i]*sint;
			R[i][q] = R[q][i] = -inMatrix[p][i]*sint + inMatrix[q][i]*cost;			 
		}
		for(i = 0;i < N; ++i)
			memcpy(inMatrix[i], R[i], sizeof(float)*N);		
		for(i = 0;i < N; ++i)
		{
			t = S[i][p];
			S[i][p] = t*cost + S[i][q]*sint;
			S[i][q] = -t*sint + S[i][q]*cost;
		}
	}

	/*Sort the eigenvector*/
	for(i = 0;i < N; ++i)
	{
		featl[i] = inMatrix[i][i];
		swap[i] = i;
	}

	for(i = 0;i < N - 1; ++i)
		for(j = i + 1;j < N; ++j)
		{
			if(featl[i] < featl[j])
			{
				t = featl[i];
				featl[i] = featl[j];
				featl[j] = t;
				sw = swap[i];
				swap[i] = swap[j];
				swap[j] = sw;
			}
		}

		for(i = 0;i < N; ++i)
			for(j = 0;j < N; ++j)
				featv[j][i] = S[j][swap[i]];

		free(swap);
		FreeArray_float(R, N, N);
		FreeArray_float(S, N, N);	
}



/* Function: fast PCA  transform, based on the book of Tsinghua, written by Bian Zhaoqi, page 224-226
 input:
	pTrainFeat   ---- Training sample(nTrainSam*ndim)
	nTrainSam:   number of sample vectors
	ndim:        the dimension of sample
 output:
	lambda       ---  EigenValue( N = min(nTrainSam,ndim) )
	trans_vector ---- EigenVector( nTrainSam*ndim )	
	total_avg    ---- Mean Feature vector
	PCA_ndim     = nTrainSam
*/
void CalPCATransParam(float *lambda, 
					  float **trans_vector,
					  float **pTrainFeat, 
					  int nTrainSam, int ndim, 
					  int PCA_ndim, float *total_avg)
{	
	float *modify_feat, *temp;
	float *c_sw;	
	int i,j,k;

	c_sw = (float *)malloc(sizeof(float)*nTrainSam*nTrainSam);		
	modify_feat = (float *)malloc(sizeof(float)*nTrainSam*ndim);
	temp = (float *)malloc(sizeof(float)*nTrainSam*ndim);

	memset(total_avg, 0, sizeof(float)*ndim);
	memset(c_sw, 0, sizeof(float)*nTrainSam*nTrainSam);

	for(i = 0;i < nTrainSam; ++i)
		for(j = 0;j < ndim; ++j)
			total_avg[j] += pTrainFeat[i][j];

	for(i = 0;i < ndim; ++i)
		total_avg[i] /= nTrainSam;

	//FILE* fp = fopen("d:\\avg.txt", "w");
	//for( k=0; k<ndim; k++)
	//	fprintf(fp, "%f ", total_avg[k]);
	//fclose(fp);

	for(i = 0;i < nTrainSam; ++i)
		for(j = 0;j < ndim; ++j)
			modify_feat[i*ndim + j] = pTrainFeat[i][j] - total_avg[j];

	MatrixTranspose(modify_feat, nTrainSam, ndim, temp);
	MatrixMulti(modify_feat, nTrainSam, ndim, temp, ndim, nTrainSam, c_sw);

	free(temp);

	float **wb;
	float **featv, *featl;

	wb = f2f(nTrainSam, nTrainSam);
	featl = (float *)malloc(sizeof(float)*nTrainSam);

	for(i = 0;i < nTrainSam; ++i)
		for(j = 0;j < nTrainSam; ++j)
			wb[i][j] = c_sw[i*nTrainSam + j];

	free(c_sw);
	featv = f2f(nTrainSam, nTrainSam);

	JacobiFeatValue(wb, featl, featv, nTrainSam);

	memset(lambda, 0, sizeof(float)*nTrainSam);
	for(i = 0;i < nTrainSam; ++i)
		lambda[i] = featl[i];

	//	float l_sum1 = 0, l_sum2 = 0;	
	//	for(i = 0;i < ndim; ++i)
	//		l_sum1 += lambda[i];
	//	for(i = 0;i < PCA_ndim; ++i)
	//		l_sum2 += lambda[i];

	for(i = 0;i < PCA_ndim; ++i)
		for(j = 0;j < ndim; ++j)
		{
			trans_vector[i][j] = 0;
			float temp_data = max(sqrt(lambda[i]), 0.000000001);
			for(k = 0;k < nTrainSam; ++k)
			{
				trans_vector[i][j] += modify_feat[k*ndim + j]*featv[k][i];
			}
			trans_vector[i][j] /= temp_data;
		}

		//	FILE *fp;
		//	fp = fopen("e:\\1.txt", "w+");
		//	fprintf(fp, "%.6f	%.6f	%.6f", l_sum1, l_sum2, l_sum2/l_sum1);
		//	fclose(fp);

		free(featl);
		free(modify_feat);	
		FreeArray_float(wb, nTrainSam, nTrainSam);
		FreeArray_float(featv, nTrainSam, nTrainSam);
}



int agmiv(double *a, int m, int n, double *b, double *x, double *aa, double eps1, double* u, double* v, int ka)
{
	int i,j;
    i=bginv(a,m,n,aa,eps1,u,v,ka);
    
	if( i<0 )
		return(-1);    
	for (i=0; i<=n-1; i++)
      { 
		x[i]=0.0;
        for(j=0; j<=m-1; j++)
          x[i]=x[i]+aa[i*m+j]*b[j];
      }
    return(1);
}


int bginv(double *a, int m, int n, double *aa, double eps1, double *u, double* v, int ka)
{ 
	int i,j,k,l,t,p,q,f;
	double wmin,wmax,thres;
	thres = 1.0e-6;
	    
    i=bmuav(a,m,n,u,v,eps1,ka);
    
	/*
	for(i=0; i<n; i++)
	{
		if(a[i*n+i]>wmax)
			wmax = a[i*n+i];
	}
	wmin = wmax*thres;
	for(i=0; i<n; i++)
	{
		if(a[i*n+i]<wmin)
			a[i*n+i] = 0;
	}
	*/

	if (i<0) 
		return(-1);
    j=n;
    if (m<n)
		j=m;
    j=j-1;
    k=0;
    while ( (k<=j) && (a[k*n+k]!=0.0) )
		k=k+1;
    k=k-1;
    for (i=0; i<=n-1; i++)
		for (j=0; j<=m-1; j++)
		{ 
			t=i*m+j; aa[t]=0.0;
			for (l=0; l<=k; l++)
			{ 
				f=l*n+i; p=j*m+l; q=l*n+l;
				aa[t]=aa[t]+v[f]*u[p]/a[q];
			}
		}
		return(1);
}
  
// a = u*a*v
int bmuav(double *a, int m, int  n, double* u, double* v, double eps1, int ka)
//int m,n,ka;
//double eps,a[],u[],v[];
{
	int i,j,k,l,it,ll,kk,ix,iy,mm,nn,iz,m1,ks;
	double d,dd,t,sm,sm1,em1,sk,ek,b,c,shh,fg[2],cs[2];
	double *s,*e,*w;
	//void ppp();
	//void sss();
	s=(double*)malloc(ka*sizeof(double));
	e=(double*)malloc(ka*sizeof(double));
	w=(double*)malloc(ka*sizeof(double));
	it=60; k=n;
	if (m-1<n) k=m-1;
	l=m;
	if (n-2<m) l=n-2;
	if (l<0) l=0;
	ll=k;
	if (l>k) ll=l;
	if (ll>=1)
	{ 
		for (kk=1; kk<=ll; kk++)
		{ 
			if (kk<=k)
			{ 
				d=0.0;
				for (i=kk; i<=m; i++)
				{ 
					ix=(i-1)*n+kk-1; 
					d=d+a[ix]*a[ix];
				}
				s[kk-1]=sqrt(d);
				if (s[kk-1]!=0.0)
				{ 
					ix=(kk-1)*n+kk-1;
					if (a[ix]!=0.0)
					{ 
						s[kk-1]=fabs(s[kk-1]);
						if (a[ix]<0.0) s[kk-1]=-s[kk-1];
					}
					for (i=kk; i<=m; i++)
					{ 
						iy=(i-1)*n+kk-1;
						a[iy]=a[iy]/s[kk-1];
					}
					a[ix]=1.0+a[ix];
				}
				s[kk-1]=-s[kk-1];
			}
			if (n>=kk+1)
			{ 
				for (j=kk+1; j<=n; j++)
				{ 
					if ((kk<=k)&&(s[kk-1]!=0.0))
					{ 
						d=0.0;
						for (i=kk; i<=m; i++)
						{ 
							ix=(i-1)*n+kk-1;
							iy=(i-1)*n+j-1;
							d=d+a[ix]*a[iy];
						}
						d=-d/a[(kk-1)*n+kk-1];
						for (i=kk; i<=m; i++)
						{ 
							ix=(i-1)*n+j-1;
							iy=(i-1)*n+kk-1;
							a[ix]=a[ix]+d*a[iy];
						}
					}
					e[j-1]=a[(kk-1)*n+j-1];
				}
			}
			if (kk<=k)
			{ 
				for (i=kk; i<=m; i++)
				{ 
					ix=(i-1)*m+kk-1; 
					iy=(i-1)*n+kk-1;
					u[ix]=a[iy];
				}
			}
			if (kk<=l)
			{ 
				d=0.0;
				for (i=kk+1; i<=n; i++)
					d=d+e[i-1]*e[i-1];
				e[kk-1]=sqrt(d);
				if (e[kk-1]!=0.0)
				{ 
					if (e[kk]!=0.0)
					{ 
						e[kk-1]=fabs(e[kk-1]);
						if (e[kk]<0.0) e[kk-1]=-e[kk-1];
					}
					for (i=kk+1; i<=n; i++)
						e[i-1]=e[i-1]/e[kk-1];
					e[kk]=1.0+e[kk];
				}
				e[kk-1]=-e[kk-1];
				if ((kk+1<=m)&&(e[kk-1]!=0.0))
				{ 
					for (i=kk+1; i<=m; i++) w[i-1]=0.0;
					for (j=kk+1; j<=n; j++)
						for (i=kk+1; i<=m; i++)
							w[i-1]=w[i-1]+e[j-1]*a[(i-1)*n+j-1];
					for (j=kk+1; j<=n; j++)
						for (i=kk+1; i<=m; i++)
						{ 
							ix=(i-1)*n+j-1;
							a[ix]=a[ix]-w[i-1]*e[j-1]/e[kk];
						}
				}
				for (i=kk+1; i<=n; i++)
					v[(i-1)*n+kk-1]=e[i-1];
			}
		}
	}
	mm=n;
	if (m+1<n) mm=m+1;
	if (k<n) s[k]=a[k*n+k];
	if (m<mm) s[mm-1]=0.0;
	if (l+1<mm) e[l]=a[l*n+mm-1];
	e[mm-1]=0.0;
	nn=m;
	if (m>n) nn=n;
	if (nn>=k+1)
	{ 
		for (j=k+1; j<=nn; j++)
		{ 
			for (i=1; i<=m; i++)
				u[(i-1)*m+j-1]=0.0;
			u[(j-1)*m+j-1]=1.0;
		}
	}
	if (k>=1)
	{ 
		for (ll=1; ll<=k; ll++)
		{ 
			kk=k-ll+1; 
			iz=(kk-1)*m+kk-1;
			if (s[kk-1]!=0.0)
			{ 
				if (nn>=kk+1)
					for (j=kk+1; j<=nn; j++)
					{ 
						d=0.0;
						for (i=kk; i<=m; i++)
						{ 
							ix=(i-1)*m+kk-1;
							iy=(i-1)*m+j-1;
							d=d+u[ix]*u[iy]/u[iz];
						}
						d=-d;
						for (i=kk; i<=m; i++)
						{ 
							ix=(i-1)*m+j-1;
							iy=(i-1)*m+kk-1;
							u[ix]=u[ix]+d*u[iy];
						}
					}
				for (i=kk; i<=m; i++)
				{
					ix=(i-1)*m+kk-1; u[ix]=-u[ix];
				}
					u[iz]=1.0+u[iz];
					if (kk-1>=1)
						for (i=1; i<=kk-1; i++)
							u[(i-1)*m+kk-1]=0.0;
			}
			else
			{ 
				for (i=1; i<=m; i++)
					u[(i-1)*m+kk-1]=0.0;
				u[(kk-1)*m+kk-1]=1.0;
			}
		}
	}
	for (ll=1; ll<=n; ll++)
	{ 
		kk=n-ll+1; 
		iz=kk*n+kk-1;
		if ((kk<=l)&&(e[kk-1]!=0.0))
		{ 
			for (j=kk+1; j<=n; j++)
			{ 
				d=0.0;
				for (i=kk+1; i<=n; i++)
				{ 
					ix=(i-1)*n+kk-1; iy=(i-1)*n+j-1;
					d=d+v[ix]*v[iy]/v[iz];
				}
				d=-d;
				for (i=kk+1; i<=n; i++)
				{ 
					ix=(i-1)*n+j-1; iy=(i-1)*n+kk-1;
					v[ix]=v[ix]+d*v[iy];
				}
			}
		}
		for (i=1; i<=n; i++)
			v[(i-1)*n+kk-1]=0.0;
		v[iz-n]=1.0;
	}
	for (i=1; i<=m; i++)
		for (j=1; j<=n; j++)
			a[(i-1)*n+j-1]=0.0;
	m1=mm; it=60;
	while (1==1)
	{ 
		if (mm==0)
		{ 
			ppp(a,e,s,v,m,n);
			free(s); 
			free(e); 
			free(w); 
			return(1);
		}
		if (it==0)
		{ 
			ppp(a,e,s,v,m,n);
			free(s); 
			free(e); 
			free(w); 
			return(-1);
		}
		kk=mm-1;
		while ((kk!=0)&&(fabs(e[kk-1])!=0.0))
		{ 
			d=fabs(s[kk-1])+fabs(s[kk]);
			dd=fabs(e[kk-1]);
			if (dd>eps1*d) kk=kk-1;
			else e[kk-1]=0.0;
		}
		if (kk==mm-1)
		{ 
			kk=kk+1;
			if (s[kk-1]<0.0)
			{ 
				s[kk-1]=-s[kk-1];
				for (i=1; i<=n; i++)
				{ 
					ix=(i-1)*n+kk-1; 
					v[ix]=-v[ix];
				}
			}
			while ((kk!=m1)&&(s[kk-1]<s[kk]))
			{ 
				d=s[kk-1]; s[kk-1]=s[kk]; s[kk]=d;
				if (kk<n)
					for (i=1; i<=n; i++)
					{ 
						ix=(i-1)*n+kk-1; iy=(i-1)*n+kk;
						d=v[ix]; v[ix]=v[iy]; v[iy]=d;
					}
				if (kk<m)
					for (i=1; i<=m; i++)
					{ 
						ix=(i-1)*m+kk-1; iy=(i-1)*m+kk;
						d=u[ix]; u[ix]=u[iy]; u[iy]=d;
					}
				kk=kk+1;
			}
			it=60;
			mm=mm-1;
		}
		else
		{ 
			ks=mm;
			while ((ks>kk)&&(fabs(s[ks-1])!=0.0))
			{ 
				d=0.0;
				if (ks!=mm) d=d+fabs(e[ks-1]);
				if (ks!=kk+1) d=d+fabs(e[ks-2]);
				dd=fabs(s[ks-1]);
				if (dd>eps1*d) ks=ks-1;
				else s[ks-1]=0.0;
			}
			if (ks==kk)
			{ 
				kk=kk+1;
				d=fabs(s[mm-1]);
				t=fabs(s[mm-2]);
				if (t>d) d=t;
				t=fabs(e[mm-2]);
				if (t>d) d=t;
				t=fabs(s[kk-1]);
				if (t>d) d=t;
				t=fabs(e[kk-1]);
				if (t>d) d=t;
				sm=s[mm-1]/d; sm1=s[mm-2]/d;
				em1=e[mm-2]/d;
				sk=s[kk-1]/d; 
				ek=e[kk-1]/d;
				b=((sm1+sm)*(sm1-sm)+em1*em1)/2.0;
				c=sm*em1; 
				c=c*c; 
				shh=0.0;
				if ((b!=0.0)||(c!=0.0))
				{ 
					shh=sqrt(b*b+c);
					if (b<0.0) shh=-shh;
					shh=c/(b+shh);
				}
				fg[0]=(sk+sm)*(sk-sm)-shh;
				fg[1]=sk*ek;
				for (i=kk; i<=mm-1; i++)
				{ 
					sss(fg,cs);
					if (i!=kk) e[i-2]=fg[0];
					fg[0]=cs[0]*s[i-1]+cs[1]*e[i-1];
					e[i-1]=cs[0]*e[i-1]-cs[1]*s[i-1];
					fg[1]=cs[1]*s[i];
					s[i]=cs[0]*s[i];
					if ((cs[0]!=1.0)||(cs[1]!=0.0))
						for (j=1; j<=n; j++)
						{ 
							ix=(j-1)*n+i-1;
							iy=(j-1)*n+i;
							d=cs[0]*v[ix]+cs[1]*v[iy];
							v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
							v[ix]=d;
						}
						sss(fg,cs);
						s[i-1]=fg[0];
						fg[0]=cs[0]*e[i-1]+cs[1]*s[i];
						s[i]=-cs[1]*e[i-1]+cs[0]*s[i];
						fg[1]=cs[1]*e[i];
						e[i]=cs[0]*e[i];
						if (i<m)
							if ((cs[0]!=1.0)||(cs[1]!=0.0))
								for (j=1; j<=m; j++)
								{ 
									ix=(j-1)*m+i-1;
									iy=(j-1)*m+i;
									d=cs[0]*u[ix]+cs[1]*u[iy];
									u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
									u[ix]=d;
								}
				}
				e[mm-2]=fg[0];
				it=it-1;
			}
			else
			{ 
				if (ks==mm)
				{ 
					kk=kk+1;
					fg[1]=e[mm-2]; 
					e[mm-2]=0.0;
					for (ll=kk; ll<=mm-1; ll++)
					{ 
						i=mm+kk-ll-1;
						fg[0]=s[i-1];
						sss(fg,cs);
						s[i-1]=fg[0];
						if (i!=kk)
						{ 
							fg[1]=-cs[1]*e[i-2];
							e[i-2]=cs[0]*e[i-2];
						}
						if ((cs[0]!=1.0)||(cs[1]!=0.0))
							for (j=1; j<=n; j++)
							{ 
								ix=(j-1)*n+i-1;
								iy=(j-1)*n+mm-1;
								d=cs[0]*v[ix]+cs[1]*v[iy];
								v[iy]=-cs[1]*v[ix]+cs[0]*v[iy];
								v[ix]=d;
							}
					}
				}
				else
				{ 
					kk=ks+1;
					fg[1]=e[kk-2];
					e[kk-2]=0.0;
					for (i=kk; i<=mm; i++)
					{ 
						fg[0]=s[i-1];
						sss(fg,cs);
						s[i-1]=fg[0];
						fg[1]=-cs[1]*e[i-1];
						e[i-1]=cs[0]*e[i-1];
						if ((cs[0]!=1.0)||(cs[1]!=0.0))
							for (j=1; j<=m; j++)
							{ 
								ix=(j-1)*m+i-1;
								iy=(j-1)*m+kk-2;
								d=cs[0]*u[ix]+cs[1]*u[iy];
								u[iy]=-cs[1]*u[ix]+cs[0]*u[iy];
								u[ix]=d;
							}
					}
				}
			}
		}
	}
	return(1);
}

static void ppp(double*a, double*e, double* s, double* v, int m, int n)
  // int m, n;
  // double a[],e[],s[],v[];
{ 
	int i,j,p,q;
	double d;
	if (m>=n) i=n;
	else i=m;
	for (j=1; j<=i-1; j++)
	{
		a[(j-1)*n+j-1]=s[j-1];
		a[(j-1)*n+j]=e[j-1];
	}
	a[(i-1)*n+i-1]=s[i-1];
	if (m<n) a[(i-1)*n+i]=e[i-1];
	for (i=1; i<=n-1; i++)
		for (j=i+1; j<=n; j++)
		{
			p=(i-1)*n+j-1; q=(j-1)*n+i-1;
			d=v[p]; v[p]=v[q]; v[q]=d;
		}
	return;
}
  
static void sss(double *fg, double* cs)
  //double cs[2],fg[2];
{ 
	double r,d;
	if ((fabs(fg[0])+fabs(fg[1]))==0.0)
	{ 
		cs[0]=1.0; 
		cs[1]=0.0; 
		d=0.0;
	}
	else 
	{ 
		d=sqrt(fg[0]*fg[0]+fg[1]*fg[1]);
		if (fabs(fg[0])>fabs(fg[1]))
		{ 
			d=fabs(d);
			if (fg[0]<0.0)
				d=-d;
		}
		if (fabs(fg[1])>=fabs(fg[0]))
		{ 
			d=fabs(d);
			if (fg[1]<0.0) d=-d;
		}
		cs[0]=fg[0]/d; 
		cs[1]=fg[1]/d;
	}
	r=1.0;
	if (fabs(fg[0])>fabs(fg[1])) 
		r=cs[1];
	else
		if (cs[0]!=0.0) r=1.0/cs[0];
	fg[0]=d; 
	fg[1]=r;
	return;
}


double* LS(double *A,double *L,int m,int n)
{
	double *At = (double*)malloc(n*m*sizeof(double));
	double *AtA= (double*)malloc(n*n*sizeof(double)); 
	double *AtL= (double*)malloc(n*sizeof(double));
	double *res= (double*)malloc(n*sizeof(double));

	transpose(A,At,m,n);
	mult(At,A,AtA,n,m,n);
	invers_matrix(AtA,n);
	mult(At,L,AtL,n,m,1);
		
	mult(AtA,AtL,res,n,n,1);
	
	free (At);
	free (AtA);
	free (AtL);

	return res;
}



//Matrix Transpose
bool MatrixTranspose(int *m1, int row1, int col1, int *m2)
{
	int i,j;
	if(m2 == NULL)
	{
		int *m3;

		m3 = (int *)malloc(sizeof(int)*row1*col1);
		for(i = 0;i < col1; ++i)
			for(j = 0;j < row1; ++j)
			{
				m3[i*col1 + j] = m1[j*col1 + i];
			}
		free(m3);
	}
	else
	{
		for(i = 0;i < col1; ++i)
			for(j = 0;j < row1; ++j)
			{
				m2[i*col1 + j] = m1[j*col1 + i];
			}
	}
	return 1;
}



bool MatrixTranspose(float *m1, int row1, int col1, float *m2)
{
	int i,j;
	if(m2 == NULL)
	{
		float *m3;

		m3 = (float *)malloc(sizeof(float)*row1*col1);
		for(i = 0;i < col1; ++i)
			for(j = 0;j < row1; ++j)
			{
				m3[i*row1 + j] = m1[j*col1 + i];
			}
		for(i = 0;i < row1; ++i)
			for(j = 0;j < col1; ++j)
				m1[i*col1 + j] = m3[j*col1 + i];
		free(m3);
	}
	else
	{
		for(i = 0;i < col1; ++i)
			for(j = 0;j < row1; ++j)
			{
				m2[i*row1 + j] = m1[j*col1 + i];
			}
	}
	return 1;
}

//Compute Inverse Matrix using complete pivot
bool MatrixInverse(float *m1, int row1, int col1)
{
	int i,j,k;
	float div,temp;
	float *out;
	int *is,*js;
	
	if(row1 != col1)
		return 0;
	
	out = (float *)malloc(sizeof(float)*row1*col1);
	is = (int *)malloc(sizeof(int)*row1);
	js = (int *)malloc(sizeof(int)*row1);
	for(i = 0;i < row1; ++i)
	{
		is[i] = i;
		js[i] = i;
	}

	// start from first column to the next 
	for(k = 0;k < row1; ++k)
	{
		div = 0;
		for(i = k;i < row1; ++i)
			for(j = k;j < row1; ++j)
			{
				if(fabs(m1[i*col1 + j]) > div)
				{
					div = fabs(m1[i*col1 + j]);
					is[k] = i;
					js[k] = j;
				}
			}
		if(fabs(div) < 1e-10)
		{
			free(out);
			free(is);
			free(js);
			return 0;
		}
		if(is[k] != k)
		{
			for(j = 0;j < row1; ++j)
			{
				temp = m1[k*col1 + j];
				m1[k*col1 + j] = m1[is[k]*col1 + j];
				m1[is[k]*col1 + j] = temp;
			}
		}
		if(js[k] != k)
		{
			for(i = 0;i < row1; ++i)
			{
				temp = m1[i*col1 + k];
				m1[i*col1 + k] = m1[i*col1 + js[k]];
				m1[i*col1 + js[k]] = temp;
			}
		}
		m1[k*col1 + k] = 1/m1[k*col1 + k];
		for(j = 0;j < row1; ++j)
		{
			if(j != k)
				m1[k*col1 + j] = m1[k*col1 + j]*m1[k*col1 + k];
		}
		for(i = 0;i < row1; ++i)
		{
			if(i != k)
			{
				for(j = 0;j < row1; ++j)
				{
					if(j != k)
						m1[i*col1 + j] -= m1[i*col1 + k]*m1[k*col1 + j];
				}
			}
		}
		for(i = 0;i < row1; ++i)
		{
			if(i != k)
				m1[i*col1 + k] = -m1[i*col1 + k]*m1[k*col1 + k];
		}							
	}

	for(k = row1 - 1;k >= 0; --k)
	{
		for(j = 0;j < row1; ++j)
			if(js[k] != k)
			{
				temp = m1[k*col1 + j];
				m1[k*col1 + j] = m1[js[k]*col1 + j];
				m1[js[k]*col1 + j] = temp;
			}
		for(i = 0;i < row1; ++i)
			if(is[k] != k)
			{
				temp = m1[i*col1 + k];
				m1[i*col1 + k] = m1[i*col1 + is[k]];
				m1[i*col1 + is[k]] = temp;
			}
	}
	free(is);
	free(js);
	free(out);
	return 1;
}



bool MatrixMulti(int *m1, int row1, int col1, int *m2, int row2, int col2,
				int *m3)
{
	int i,j,k;

	if(col1 != row2)
		return 0;
	for(i = 0;i < row1; ++i)
		for(j = 0;j < col2; ++j)
		{
			int sum = 0;
			for(k = 0;k < col1; ++k)
				sum += m1[i*col1 + k]*m2[k*col2 + j];
			m3[i*col2 + j] = sum;
		}
	
	return 1;
}


bool MatrixMulti(float *m1, int row1, int col1, float *m2, int row2, int col2,
				float *m3)
{
	int i,j,k;

	if(col1 != row2)
		return 0;
	for(i = 0;i < row1; ++i)
		for(j = 0;j < col2; ++j)
		{
			float sum = 0;
			for(k = 0;k < col1; ++k)
				sum += m1[i*col1 + k]*m2[k*col2 + j];
			m3[i*col2 + j] = sum;
		}
	
	return 1;
}

bool MatrixMulti(float *m1, int row1, int col1, int *m2, int row2, int col2,
				float *m3)
{
	int i,j,k;

	if(col1 != row2)
		return 0;
	for(i = 0;i < row1; ++i)
		for(j = 0;j < col2; ++j)
		{
			float sum = 0;
			for(k = 0;k < col1; ++k)
				sum += m1[i*col1 + k]*m2[k*col2 + j];
			m3[i*col2 + j] = sum;
		}
	
	return 1;
}

bool MatrixMulti(int *m1, int row1, int col1, float *m2, int row2, int col2,
				float *m3)
{
	int i,j,k;

	if(col1 != row2)
		return 0;
	for(i = 0;i < row1; ++i)
		for(j = 0;j < col2; ++j)
		{
			float sum = 0;
			for(k = 0;k < col1; ++k)
				sum += m1[i*col1 + k]*m2[k*col2 + j];
			m3[i*col2 + j] = sum;
		}
	
	return 1;
}