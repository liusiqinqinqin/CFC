/***********************************************************************************
Code for difference image generation methods, including log-ratio, mean-ratio, and neighborhood-ratio methods.

3 functions:
log_ratio: code for log-ratio method;
NR: code for neighborhood-ratio method;
mean_ratio: code for mean-ratio method.

input:
x1: a matrix, the first image with the data type of double;
x2: a matrix, the second image with the data type of double;
width: the width of the two image;
height: the height of the two image;
wn: size of the neighborhood (wn*wn) in neighborhood-ratio and mean-ratio methods.

output:
the return of the function: a matrix, the difference image with the data type of double.

***********************************************************************************/
#include<stdio.h>
#include<iostream>
#include<malloc.h>
#include<math.h>
#include<highgui.h>
#include<time.h>





#define e 2.7182882846
#define pi 3.141592657


void cal_gm(double **f, double diat, int n)                         //计算f(x,y),高斯滤波模板
{
	int i, j;
	double s;
	s = 0;
	for (j = 0; j<n / 2 + 1; j++)
	{
		for (i = 0; i<n / 2 + 1; i++)
		{
			f[n / 2 - i][n / 2 - j] = pow(e, -(i*i + j*j) / (2 * diat*diat));
			s = s + f[n / 2 - i][n / 2 - j];
		}
	}
	for (i = 0; i<n / 2 + 1; i++)
	{
		for (j = n / 2 + 1; j<n; j++)
		{
			f[i][j] = f[i][n - 1 - j];
			s = s + f[i][j];
		}
	}
	for (i = n / 2 + 1; i<n; i++)
	{
		for (j = 0; j<n / 2 + 1; j++)
		{
			f[i][j] = f[n - 1 - i][j];
			s = s + f[i][j];
		}
	}
	for (i = n / 2 + 1; i<n; i++)
	{
		for (j = n / 2 + 1; j<n; j++)
		{
			f[i][j] = f[n - 1 - i][n - 1 - j];
			s = s + f[i][j];
		}
	}


	for (i = 0; i<n; i++)
	{
		for (j = 0; j<n; j++)
		{
			f[i][j] = f[i][j] / s;
		}
	}
}




void g_filter(double **x, double**y, int width, int height, int n)
{
	int i, j, k, l;
	double **gm;
	gm = (double**)malloc(n * sizeof(double*));
	for (i = 0; i<n; i++)
	{
		gm[i] = (double*)malloc(n * sizeof(double));
	}

	cal_gm(gm, 1.3, n);


	double s;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			s = 0;
			for (k = 0; k<n; k++)
			{
				for (l = 0; l<n; l++)
				{
					if (i - n / 2 + k >= 0 && i - n / 2 + k<width&&j - n / 2 + l >= 0 && j - n / 2 + l<height)
					{
						s += gm[k][l] * x[i - n / 2 + k][j - n / 2 + l];
					}
				}
			}
			y[i][j] = s;
		}
	}

	for (i = 0; i<n; i++)
	{
		free(gm[i]);
	}
	free(gm);
}


void cal_lp(double **x, double ***lp, int width, int height, int l)
{
	int i, j, k;
	double ***g;
	g = (double***)malloc(l * sizeof(double**));
	for (k = 0; k<l; k++)
	{
		g[k] = (double**)malloc(width * sizeof(double*));
		for (i = 0; i<width; i++)
		{
			g[k][i] = (double*)malloc(height * sizeof(double));
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			g[0][i][j] = x[i][j];
		}
	}

	for (k = 1; k<l; k++)
	{
		g_filter(g[k - 1], g[k], width, height, 5);
	}


	for (k = 0; k<l - 1; k++)
	{
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				lp[k][i][j] = g[k][i][j] - g[k + 1][i][j];
			}
		}
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			lp[l - 1][i][j] = g[l - 1][i][j];
		}
	}

	for (k = 0; k<l; k++)
	{
		for (i = 0; i<width; i++)
		{
			free(g[k][i]);
		}
		free(g[k]);
	}
	free(g);
}



void reconstruct_lp(double ***lp, double **y, int width, int height, int l)
{
	int i, j, k;
	for (k = l - 2; k >= 0; k--)
	{
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				lp[k][i][j] = lp[k][i][j] + lp[k + 1][i][j];
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = lp[0][i][j];
		}
	}
}







void filter(double **x, double **y, double **model, int width, int height, int n)
{
	int i, j, k, l;
	double s;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			s = 0;
			for (k = 0; k<n; k++)
			{
				for (l = 0; l<n; l++)
				{
					if (i - n / 2 + k >= 0 && i - n / 2 + k<width&&j - n / 2 + l >= 0 && j - n / 2 + l<height)
					{
						s += model[k][l] * x[i - n / 2 + k][j - n / 2 + l];//printf("%f\n",model[k][l]*x[i-n/2+k][j-n/2+l]);
					}
				}
			}
			y[i][j] = s;
		}
	}
}



void direction4_filter(double **x, double ***y, int width, int height)
{
	int i, j, k;
	double ***mod;

	int wn = 11;

	mod = (double***)malloc(4 * sizeof(double**));
	for (i = 0; i<4; i++)
	{
		mod[i] = (double**)malloc(wn * sizeof(double*));
		for (j = 0; j<wn; j++)
		{
			mod[i][j] = (double*)malloc(wn * sizeof(double));
		}
	}

	FILE *data;
	fopen_s(&data,"data4.txt", "r+");
	for (k = 0; k<4; k++)
	{
		for (i = 0; i<wn; i++)
		{
			for (j = 0; j<wn; j++)
			{
				fscanf_s(data, "%lf", &mod[k][i][j]);
			}
		}
	}
	fclose(data);

	for (k = 0; k<4; k++)
	{
		filter(x, y[k], mod[k], width, height, wn);
	}
	for (k = 0; k<4; k++)
	{
		for (i = 0; i<wn; i++)
		{
			free(mod[k][i]);
		}
		free(mod[k]);
	}
	free(mod);
}


void direction8_filter(double**x, double ***y, int width, int height)
{
	int i, j, k;
	double ***mod;

	int wn = 11;

	mod = (double***)malloc(8 * sizeof(double**));
	for (k = 0; k<8; k++)
	{
		mod[k] = (double**)malloc(wn * sizeof(double*));
		for (i = 0; i<wn; i++)
		{
			mod[k][i] = (double*)malloc(wn * sizeof(double));
		}
	}

	FILE *data;
	fopen_s(&data,"data8.txt", "r+");
	for (k = 0; k<8; k++)
	{
		for (i = 0; i<wn; i++)
		{
			for (j = 0; j<wn; j++)
			{
				fscanf_s(data, "%lf", &mod[k][i][j]);
			}
		}
	}
	fclose(data);


	for (k = 0; k<8; k++)
	{
		filter(x, y[k], mod[k], width, height, wn);
	}


	for (k = 0; k<8; k++)
	{
		for (i = 0; i<wn; i++)
		{
			free(mod[k][i]);
		}
		free(mod[k]);
	}
	free(mod);
}



void reconstruct_direction4(double ***x, double**y, int width, int height)
{
	int i, j, k;
	/*double ***mod;

	int wn=11;

	mod=(double***)malloc(4*sizeof(double**));
	for(k=0;k<4;k++)
	{
	mod[k]=(double**)malloc(wn*sizeof(double*));
	for(i=0;i<wn;i++)
	{
	mod[k][i]=(double*)malloc(wn*sizeof(double));
	}
	}

	FILE *data;
	data=fopen("data4.txt","r+");
	for(k=0;k<4;k++)
	{
	for(i=0;i<wn;i++)
	{
	for(j=0;j<wn;j++)
	{
	fscanf_s(data,"%lf",&mod[k][i][j]);
	}
	}
	}
	fclose(data);

	double ***mid;
	mid=(double***)malloc(4*sizeof(double**));
	for(k=0;k<4;k++)
	{
	mid[k]=(double**)malloc(width*sizeof(double*));
	for(i=0;i<width;i++)
	{
	mid[k][i]=(double*)malloc(height*sizeof(double));
	}
	}

	for(k=0;k<4;k++)
	{
	filter(x[k],mid[k],mod[k],width,height,wn);
	}*/

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 0;
			for (k = 0; k<4; k++)
			{
				//y[i][j]+=mid[k][i][j];
				y[i][j] += x[k][i][j];
			}
		}
	}

	/*for(k=0;k<4;k++)
	{
	for(i=0;i<wn;i++)
	{
	free(mod[k][i]);
	}
	free(mod[k]);
	for(i=0;i<width;i++)
	{
	free(mid[k][i]);
	}
	free(mid[k]);
	}
	free(mod);
	free(mid);*/
}


void reconstruct_direction8(double ***x, double**y, int width, int height)
{
	int i, j, k;
	/*double ***mod;

	int wn=11;

	mod=(double***)malloc(8*sizeof(double**));
	for(k=0;k<8;k++)
	{
	mod[k]=(double**)malloc(wn*sizeof(double*));
	for(i=0;i<wn;i++)
	{
	mod[k][i]=(double*)malloc(wn*sizeof(double));
	}
	}

	FILE *data;
	data=fopen("data8.txt","r+");
	for(k=0;k<8;k++)
	{
	for(i=0;i<wn;i++)
	{
	for(j=0;j<wn;j++)
	{
	fscanf_s(data,"%lf",&mod[k][i][j]);
	}
	}
	}
	fclose(data);

	double ***mid;
	mid=(double***)malloc(8*sizeof(double**));
	for(k=0;k<8;k++)
	{
	mid[k]=(double**)malloc(width*sizeof(double*));
	for(i=0;i<width;i++)
	{
	mid[k][i]=(double*)malloc(height*sizeof(double));
	}
	}

	for(k=0;k<8;k++)
	{
	filter(x[k],mid[k],mod[k],width,height,wn);
	}*/

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 0;
			for (k = 0; k<8; k++)
			{
				//y[i][j]+=mid[k][i][j];
				y[i][j] += x[k][i][j];
			}
		}
	}

	/*for(k=0;k<8;k++)
	{
	for(i=0;i<wn;i++)
	{
	free(mod[k][i]);
	}
	free(mod[k]);
	for(i=0;i<width;i++)
	{
	free(mid[k][i]);
	}
	free(mid[k]);
	}
	free(mod);
	free(mid);*/
}



void nsct(double **x, double **low, double ***high4, double ***high8, int width, int height)
{
	int i, j;
	double ***lp;
	lp = (double***)malloc(3 * sizeof(double**));
	for (i = 0; i<3; i++)
	{
		lp[i] = (double**)malloc(width * sizeof(double*));
		for (j = 0; j<width; j++)
		{
			lp[i][j] = (double*)malloc(height * sizeof(double));
		}
	}

	cal_lp(x, lp, width, height, 3);

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			low[i][j] = lp[2][i][j];
		}
	}

	direction4_filter(lp[1], high4, width, height);
	direction8_filter(lp[0], high8, width, height);

	for (i = 0; i<3; i++)
	{
		for (j = 0; j<width; j++)
		{
			free(lp[i][j]);
		}
		free(lp[i]);
	}
	free(lp);
}


void insct(double**low, double ***high4, double ***high8, double**y, int width, int height)
{
	int i, j;

	double ***lp;
	lp = (double***)malloc(3 * sizeof(double**));
	for (i = 0; i<3; i++)
	{
		lp[i] = (double**)malloc(width * sizeof(double*));
		for (j = 0; j<width; j++)
		{
			lp[i][j] = (double*)malloc(height * sizeof(double));
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			lp[2][i][j] = low[i][j];
		}
	}

	reconstruct_direction4(high4, lp[1], width, height);
	reconstruct_direction8(high8, lp[0], width, height);

	reconstruct_lp(lp, y, width, height, 3);

	for (i = 0; i<3; i++)
	{
		for (j = 0; j<width; j++)
		{
			free(lp[i][j]);
		}
		free(lp[i]);
	}
	free(lp);
}









/*double ek(double**im,int x,int y,long width,long height,int n)
{
int i,j;
double s;
s=0;
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
if(x-n/2+i>=0&&x-n/2+i<width&&y-n/2+j>=0&&y-n/2+j<height)
{
s+=(im[x-n/2+i][y-n/2+j]-im[x][y])*(im[x-n/2+i][y-n/2+j]-im[x][y]);
}
}
}
return s/(im[x][y]*im[x][y]);
}


double locale(double **im,int x,int y,long width,long height,double **w,int n)
{
int i,j;
double s;
s=0;
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
if(x-n/2+i>=0&&x-n/2+i<width&&y-n/2+j>=0&&y-n/2+j<height)
{
s+=w[i][j]*(im[x-n/2+i][y-n/2+j]-im[x][y])*(im[x-n/2+i][y-n/2+j]-im[x][y]);
}
}
}
return s;
}



double ab(double x)
{
if(x>0)
{
return x;
}
else
{
return -x;
}
}

double matmatrix(double **im1,double **im2,int x,int y,long width,long height,double **w,int n)
{
int i,j;
double s;
s=0;
for(i=0;i<n;i++)
{
for(j=0;j<n;j++)
{
if(x-n/2+i>=0&&x-n/2+i<width&&y-n/2+j>=0&&y-n/2+j<height)
{
s+=w[i][j]*ab(im1[x-n/2+i][y-n/2+j]-im1[x][y])*ab(im2[x-n/2+i][y-n/2+j]-im2[x][y]);
}
}
}
return 2*s/(locale(im1,x,y,width,height,w,n)+locale(im2,x,y,width,height,w,n));
}



double sse(double **x,int i,int j,long width,long height,int wn)
{
int k,l;
double u,diat,num;
num=0;u=0;diat=0;
for(k=0;k<wn;k++)
{
for(l=0;l<wn;l++)
{
if(i-wn/2+k>=0&&i-wn/2+k<width&&j-wn/2+l>=0&&j-wn/2+l<height)
{
u+=x[i-wn/2+k][j-wn/2+l];
num++;
}
}
}
u=u/num;
for(k=0;k<wn;k++)
{
for(l=0;l<wn;l++)
{
if(i-wn/2+k>=0&&i-wn/2+k<width&&j-wn/2+l>=0&&j-wn/2+l<height)
{
diat+=(x[i-wn/2+k][j-wn/2+l]-u)*(x[i-wn/2+k][j-wn/2+l]-u);
}
}
}
diat=sqrt(diat/num);
return diat;
}
*/

/*void fusion_rule_low(double **x1,double **x2,double **y,long width,long height)
{
int i,j;
double q=0.49;

for(i=0;i<width;i++)
{
for(j=0;j<height;j++)
{
/*if(x1[i][j]==0||x2[i][j]==0)
{
y[i][j]=0.5*x1[i][j]+0.5*x2[i][j];
}
else if(ek(x1,i,j,width,height,5)<ek(x2,i,j,width,height,5))
{
y[i][j]=x2[i][j]*(1-q)+x1[i][j]*q;
}
else
{
y[i][j]=x1[i][j]*(1-q)+x2[i][j]*q;
}*/
/*if(sse(x1,i,j,width,height,5)<sse(x2,i,j,width,height,5))
{
y[i][j]=x1[i][j]*(1-q)+x2[i][j]*q;
}
else
{
y[i][j]=x2[i][j]*(1-q)+x1[i][j]*q;
}
}
}
}



void fusion_rule_high(double **x1,double **x2,double **y,long width,long height)
{
int i,j;
double **w;
int nw=3;
w=(double**)malloc(nw*sizeof(double*));
for(i=0;i<nw;i++)
{
w[i]=(double*)malloc(nw*sizeof(double));
}
cal_gm(w,15,nw);

double t=0.57;
double q=0;

for(i=0;i<width;i++)
{
for(j=0;j<height;j++)
{
if(matmatrix(x1,x2,i,j,width,height,w,nw)<t)
{
if(x1[i][j]>x2[i][j])
//if(locale(x1,i,j,width,height,w,nw)>locale(x2,i,j,width,height,w,nw))
//if(sse(x1,i,j,width,height,5)<sse(x2,i,j,width,height,5))
{
y[i][j]=x1[i][j]*(1-q)+x2[i][j]*q;
}
else
{
y[i][j]=x2[i][j]*(1-q)+x1[i][j]*q;
}
}
else
{
if(x1[i][j]<x2[i][j])
//if(locale(x1,i,j,width,height,w,nw)<locale(x2,i,j,width,height,w,nw))
//if(sse(x1,i,j,width,height,5)>sse(x2,i,j,width,height,5))
{
y[i][j]=x1[i][j]*(1-q)+x2[i][j]*q;
}
else
{
y[i][j]=x2[i][j]*(1-q)+x1[i][j]*q;
}
}
/*if(sse(x1,i,j,width,height,5)<sse(x2,i,j,width,height,5))
{
y[i][j]=x1[i][j]*(1-q)+x2[i][j]*q;
}
else
{
y[i][j]=x2[i][j]*(1-q)+x1[i][j]*q;
}*/
/*}
}

for(i=0;i<nw;i++)
{
free(w[i]);
}
free(w);
}
*/


void fusion_rule_low(double **x1, double **x2, double **y, int width, int height)
{
	int i, j;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			//y[i][j]=0.2*x2[i][j]+0.8*x1[i][j];
			if (x1[i][j]<x2[i][j])
			{
				y[i][j] = x1[i][j];
			}
			else
			{
				y[i][j] = x2[i][j];
			}
		}
	}
}

void fusion_rule_high(double **x1, double **x2, double **y, int width, int height)
{
	int i, j, k, l;
	double s1, s2;
	int wn = 3;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			s1 = 0; s2 = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						s1 += x1[i - wn / 2 + k][j - wn / 2 + l] * x1[i - wn / 2 + k][j - wn / 2 + l];
						s2 += x2[i - wn / 2 + k][j - wn / 2 + l] * x2[i - wn / 2 + k][j - wn / 2 + l];
					}
				}
			}
			if (s1<s2)
			{
				y[i][j] = x1[i][j];
			}
			else
			{
				y[i][j] = x2[i][j];
			}
		}
	}
}




void nsct_fusion(double **x1, double **x2, double **y, int width, int height)
{
	int i, k;
	double **low1, **low2, **low;
	double ***high41, ***high42, ***high4;
	double ***high81, ***high82, ***high8;

	low1 = (double**)malloc(width * sizeof(double*));
	low2 = (double**)malloc(width * sizeof(double*));
	low = (double**)malloc(width * sizeof(double*));
	for (i = 0; i<width; i++)
	{
		low1[i] = (double*)malloc(height * sizeof(double));
		low2[i] = (double*)malloc(height * sizeof(double));
		low[i] = (double*)malloc(height * sizeof(double));
	}

	high41 = (double***)malloc(4 * sizeof(double**));
	high42 = (double***)malloc(4 * sizeof(double**));
	high4 = (double***)malloc(4 * sizeof(double**));
	for (k = 0; k<4; k++)
	{
		high41[k] = (double**)malloc(width * sizeof(double*));
		high42[k] = (double**)malloc(width * sizeof(double*));
		high4[k] = (double**)malloc(width * sizeof(double*));
		for (i = 0; i<width; i++)
		{
			high41[k][i] = (double*)malloc(height * sizeof(double));
			high42[k][i] = (double*)malloc(height * sizeof(double));
			high4[k][i] = (double*)malloc(height * sizeof(double));
		}
	}


	high81 = (double***)malloc(8 * sizeof(double**));
	high82 = (double***)malloc(8 * sizeof(double**));
	high8 = (double***)malloc(8 * sizeof(double**));
	for (k = 0; k<8; k++)
	{
		high81[k] = (double**)malloc(width * sizeof(double*));
		high82[k] = (double**)malloc(width * sizeof(double*));
		high8[k] = (double**)malloc(width * sizeof(double*));
		for (i = 0; i<width; i++)
		{
			high81[k][i] = (double*)malloc(height * sizeof(double));
			high82[k][i] = (double*)malloc(height * sizeof(double));
			high8[k][i] = (double*)malloc(height * sizeof(double));
		}
	}


	nsct(x1, low1, high41, high81, width, height);
	printf("decomposing image 1 complete-------------\n");

	nsct(x2, low2, high42, high82, width, height);
	printf("decomposing image 2 complete-------------\n");

	fusion_rule_low(low1, low2, low, width, height);
	printf("fusing lowpass image complete------------\n");

	for (k = 0; k<4; k++)
	{
		fusion_rule_high(high41[k], high42[k], high4[k], width, height);
	}

	for (k = 0; k<8; k++)
	{
		fusion_rule_high(high81[k], high82[k], high8[k], width, height);
	}
	printf("fusing highpass image complete-----------\n");

	insct(low, high4, high8, y, width, height);
	printf("reconstruction complete------------------\n");




	for (k = 0; k<4; k++)
	{
		for (i = 0; i<width; i++)
		{
			free(high41[k][i]);
			free(high42[k][i]);
			free(high4[k][i]);
		}
		free(high41[k]);
		free(high42[k]);
		free(high4[k]);
	}
	free(high41);
	free(high42);
	free(high4);

	for (k = 0; k<8; k++)
	{
		for (i = 0; i<width; i++)
		{
			free(high81[k][i]);
			free(high82[k][i]);
			free(high8[k][i]);
		}
		free(high81[k]);
		free(high82[k]);
		free(high8[k]);
	}
	free(high81);
	free(high82);
	free(high8);

	for (i = 0; i<width; i++)
	{
		free(low1[i]);
		free(low2[i]);
		free(low[i]);
	}
	free(low1);
	free(low2);
	free(low);
}







