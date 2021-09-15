/******************************************************************************************************************************************************
code for wavelet fusion method which is used to generate the difference image.

1main function:
wavelet_fusion: code for generating the difference image by wavelet_fusion method.

input:
x1: a matrix, the first image with the data type of double;
x2: a matrix, the second image with the data type of double;
width: the width of the two image;
height: the height of the two image;
n: the number of decomposition levels in wavelet decomposition.

output:
the return of the function: a matrix, the difference image with the data type of double.

The rest functions and variables are intermediate functions which are steps of wavelet_fusion and parameters of filters.
******************************************************************************************************************************************************/


#include"SAR.h"


#define pi 3.141592657




double h[10] = { 0.0033, -0.0126, -0.0062, 0.0776, -0.0322, -0.2423, 0.1384, 0.7243, 0.6038, 0.1601 };
double g[10] = { -0.1601, 0.6038, -0.7243, 0.1384, 0.2423, -0.0322, -0.0776, -0.0062, 0.0126, 0.0033 };
int N = 10;


void row(double *x, double *y1, double *y2, int length, int n)      //length x;n g
{
	int i, j;
	double s1, s2;
	for (i = 0; i<length; i++)
	{
		s1 = 0; s2 = 0;
		for (j = 0; j<n; j++)
		{
			s1 += h[j] * x[(i + j) % length];
			s2 += g[j] * x[(i + j) % length];
		}
		y1[i] = s1;
		y2[i] = s2;
	}
}

void wavelet_tra(double **x, double **ll, double **hl, double **lh, double **hh, long width, long height)
{
	int i, j;
	double *mid, *y1, *y2;
	mid = (double*)malloc(width*sizeof(double));
	y1 = (double*)malloc(width*sizeof(double));
	y2 = (double*)malloc(width*sizeof(double));
	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			mid[j] = x[j][i];
		}
		row(mid, y1, y2, width, N);
		for (j = 0; j<width; j++)
		{
			ll[j][i] = y1[j];
			hl[j][i] = y2[j];
		}
	}
	free(mid);
	free(y1);
	free(y2);
	mid = (double*)malloc(height*sizeof(double));
	y1 = (double*)malloc(height*sizeof(double));
	y2 = (double*)malloc(height*sizeof(double));
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid[j] = ll[i][j];
		}
		row(mid, y1, y2, height, N);
		for (j = 0; j<height; j++)
		{
			ll[i][j] = y1[j];
			lh[i][j] = y2[j];
		}
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid[j] = hl[i][j];
		}
		row(mid, y1, y2, height, N);
		for (j = 0; j<height; j++)
		{
			hl[i][j] = y1[j];
			hh[i][j] = y2[j];
		}
	}
	free(mid);
	free(y1);
	free(y2);
}


void irow(double *low, double *high, double *y, int length, int n)
{
	int i, j;
	double s1, s2;
	for (i = 0; i<length; i++)
	{
		s1 = 0; s2 = 0;
		for (j = 0; j<n; j++)
		{
			s1 += h[j] * low[length - 1 - ((length - 1 - i + j) % length)];
			s2 += g[j] * high[length - 1 - ((length - 1 - i + j) % length)];
		}
		y[i] = s1 + s2;
	}
}




void iwavelet_tra(double **ll, double **hl, double **lh, double **hh, double **y, long width, long height)
{
	int i, j;
	double *mid, *low, *high;
	mid = (double*)malloc(height*sizeof(double));
	high = (double*)malloc(height*sizeof(double));
	low = (double*)malloc(height*sizeof(double));
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			low[j] = ll[i][j];
			high[j] = lh[i][j];
		}
		irow(low, high, mid, height, N);
		for (j = 0; j<height; j++)
		{
			ll[i][j] = mid[j];
		}
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			low[j] = hl[i][j];
			high[j] = hh[i][j];
		}
		irow(low, high, mid, height, N);
		for (j = 0; j<height; j++)
		{
			hl[i][j] = mid[j];
		}
	}
	free(mid);
	free(low);
	free(high);
	mid = (double*)malloc(width*sizeof(double));
	low = (double*)malloc(width*sizeof(double));
	high = (double*)malloc(width*sizeof(double));
	for (i = 0; i<height; i++)
	{
		for (j = 0; j<width; j++)
		{
			low[j] = ll[j][i];
			high[j] = hl[j][i];
		}
		irow(low, high, mid, width, N);
		for (j = 0; j<width; j++)
		{
			y[j][i] = mid[j];
		}
	}
	free(low);
	free(high);
	free(mid);
}






void dwt(double **x, double ***y, long width, long height, int n)
{
	int i, j, k;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[3 * n][i][j] = x[i][j];
		}
	}
	for (k = n - 1; k >= 0; k--)
	{
		wavelet_tra(y[3 * (k + 1)], y[3 * k], y[3 * k + 1], y[3 * k + 2], y[3 * k + 3], width, height);
	}
}

void idwt(double ***x, double **y, long width, long height, int n)
{
	int i, j, k;
	for (k = 0; k<n; k++)
	{
		iwavelet_tra(x[k * 3], x[k * 3 + 1], x[k * 3 + 2], x[k * 3 + 3], x[k * 3 + 3], width, height);
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = x[3 * n][i][j];
		}
	}
}




double** wavelet_fusion(double **x1, double **x2, long width, long height, int n)
{
	double ***wavelet1, ***wavelet2, **y;
	int i, j, k;
	y = new double*[width];
	for (i = 0; i < width; i++)
	{
		y[i] = new double[height];
	}

	double**lr, **mr;


	lr = log_ratio(x1, x2, width, height);
	mr = mean_ratio(x1, x2, width, height, 3);

	wavelet1 = (double***)malloc((3 * n + 1)*sizeof(double**));
	wavelet2 = (double***)malloc((3 * n + 1)*sizeof(double**));
	for (k = 0; k<3 * n + 1; k++)
	{
		wavelet1[k] = (double**)malloc(width*sizeof(double*));
		wavelet2[k] = (double**)malloc(width*sizeof(double*));
		for (i = 0; i<width; i++)
		{
			wavelet1[k][i] = (double*)malloc(height*sizeof(double));
			wavelet2[k][i] = (double*)malloc(height*sizeof(double));
		}
	}
	dwt(lr, wavelet1, width, height, n);
	dwt(mr, wavelet2, width, height, n);
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			wavelet1[0][i][j] = wavelet2[0][i][j];
		}
	}
	idwt(wavelet1, y, width, height, n);

	for (i = 0; i < 3 * n + 1; i++)
	{
		for (j = 0; j < width; j++)
		{
			delete[]wavelet1[i][j];
			delete[]wavelet2[i][j];
		}
		delete[]wavelet1[i];
		delete[]wavelet2[i];
	}
	delete[]wavelet1;
	delete[]wavelet2;


	double maximum, minimum;
	maximum = 0; minimum = 10000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (y[i][j]>maximum)
			{
				maximum = y[i][j];
			}
			if (y[i][j]<minimum)
			{
				minimum = y[i][j];
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 255.0*(y[i][j] - minimum) / (maximum - minimum);
		}
	}

	for (i = 0; i < width; i++)
	{
		delete[]lr[i];
		delete[]mr[i];
	}
	delete[]lr;
	delete[]mr;

	return y;
}
