/************************************************************************************************************************************
code for recognizing changed regions from the difference image, including FLICM and RFLICM.

2 main fuctions:
flicm: code for FLICM method;
rflicm: code for RFLICM;

input:
x: a matrix, the difference image with the data type of double;
width: the width of the difference image;
height: the height of the difference image;
m: the clustering parameter, is set to 2;
wn: the size of the neighborhood (wn*wn).

output:
the return of the fuction: a matrix, the change detection map, with white pixels the changed pixel.




*************************************************************************************************************************************/



#include"SAR.h"

void cal_v(double **x, int width, int height, double **u1, double **u2, double &v1, double &v2, int m)
{
	double s1, s11, s2, s21, mid1, mid2;
	s1 = 0; s11 = 0; s2 = 0; s21 = 0;
	int i, j;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid1 = pow(u1[i][j], m);
			mid2 = pow(u2[i][j], m);
			s1 += mid1*x[i][j] / 100;
			s2 += mid2*x[i][j] / 100;
			s11 += mid1;
			s21 += mid2;
		}
	}
	s11 = s11 / 100;
	s21 = s21 / 100;

	v1 = s1 / s11;
	v2 = s2 / s21;
}


void cal_u(double **x, int width, int height, double v1, double v2, double **u1, double **u2, double **wij, int wn, int m)
{
	double **g1, **g2;
	int i, j;
	g1 = new double*[width];
	g2 = new double*[width];
	for (i = 0; i<width; i++)
	{
		g1[i] = new double[height];
		g2[i] = new double[height];
	}

	int k, l, px, py;
	double num;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			g1[i][j] = 0;
			g2[i][j] = 0;
			num = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k; py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						g1[i][j] += wij[k][l] * pow(1 - u1[px][py], m)*(x[px][py] - v1)*(x[px][py] - v1);
						g2[i][j] += wij[k][l] * pow(1 - u2[px][py], m)*(x[px][py] - v2)*(x[px][py] - v2);
						num++;
					}
				}
			}
			//printf("%f  %f\n",g1[i][j],g2[i][j]);
			//g1[i][j]=g1[i][j]/num;
			//g2[i][j]=g2[i][j]/num;
		}
	}

	double mid;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid = pow(((x[i][j] - v1)*(x[i][j] - v1) + g1[i][j]) / ((x[i][j] - v2)*(x[i][j] - v2) + g2[i][j]), 1 / (m - 1)) + 1;
			u1[i][j] = 1 / mid;
			mid = pow(((x[i][j] - v2)*(x[i][j] - v2) + g2[i][j]) / ((x[i][j] - v1)*(x[i][j] - v1) + g1[i][j]), 1 / (m - 1)) + 1;
			u2[i][j] = 1 / mid;
		}
	}

	for (i = 0; i<width; i++)
	{
		delete[]g1[i];
		delete[]g2[i];
	}
	delete[]g1;
	delete[]g2;
}

void initial_u(double **x, int width, int height, double **u1, double **u2, int thre)
{
	srand((unsigned int)time(NULL));
	int i, j;
	double mid;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid = (double)(rand() % 10000) / 10000.0;
			if (x[i][j]>thre)
			{
				if (mid<0.5)
				{
					u1[i][j] = 1 - mid;
					u2[i][j] = mid;
				}
				else
				{
					u1[i][j] = mid;
					u2[i][j] = 1 - mid;
				}
			}
			else
			{
				if (mid<0.5)
				{
					u1[i][j] = mid;
					u2[i][j] = 1 - mid;
				}
				else
				{
					u1[i][j] = 1 - mid;
					u2[i][j] = mid;
				}
			}
		}
	}
}
double** flicm(double **x, int width, int height, int m, int wn)
{
	int i, j;
	double **wij;
	wij = new double*[wn];
	for (i = 0; i<wn; i++)
	{
		wij[i] = new double[wn];
	}

	double dis;
	for (i = 0; i<wn; i++)
	{
		for (j = 0; j<wn; j++)
		{
			dis = sqrt((double)((i - (int)(wn / 2))*(i - (int)(wn / 2))) + (double)((j - (int)(wn / 2))*(j - (int)(wn / 2))));
			wij[i][j] = 1 / (dis + 1);
		}
	}
	wij[wn / 2][wn / 2] = 0;

	double **u1, **u2;
	u1 = new double*[width];
	u2 = new double*[width];
	for (i = 0; i<width; i++)
	{
		u1[i] = new double[height];
		u2[i] = new double[height];
	}

	initial_u(x, width, height, u1, u2, 150);

	double v1, v2, midv1, midv2;
	midv1 = 0; midv2 = 255;
	cal_v(x, width, height, u1, u2, v1, v2, m);
	while ((v1 - midv1)*(v1 - midv1) + (v2 - midv2)*(v2 - midv2)>0.001)
	{
		cal_u(x, width, height, v1, v2, u1, u2, wij, wn, m);
		midv1 = v1; midv2 = v2;
		cal_v(x, width, height, u1, u2, v1, v2, m);
	}

	cal_u(x, width, height, v1, v2, u1, u2, wij, wn, m);

	double **y;
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (u1[i][j] >= u2[i][j])
			{
				y[i][j] = 0;
			}
			else
			{
				y[i][j] = 255;
			}
		}
	}

	for (i = 0; i<wn; i++)
	{
		delete[]wij[i];
	}
	delete[]wij;
	for (i = 0; i<width; i++)
	{
		delete[]u1[i];
		delete[]u2[i];
	}
	delete[]u1;
	delete[]u2;
	return y;
}

double max(double x1, double x2)
{
	if (x1>x2)
	{
		return x1;
	}
	else
	{
		return x2;
	}
}
void local1(double **x, double **&y, int width, int height, int wn, double lmds, double lmdg)
{
	int i, j, k, l;
	double diat, s1, s2;
	int N;
	double *S;
	S = new double[wn*wn];

	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			N = 0; diat = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						diat += (x[i][j] - x[i - wn / 2 + k][j - wn / 2 + l])*(x[i][j] - x[i - wn / 2 + k][j - wn / 2 + l]);
						N++;
					}
				}
			}
			diat = diat / (double)(N - 1);

			N = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						if (i == i - wn / 2 + k&&j == j - wn / 2 + l)
						{
							S[N] = 0;
							N++;
						}
						else
						{
							S[N] = exp(-max((double)abs(k - wn / 2), (double)abs(l - wn / 2)) / lmds - (x[i][j] - x[i - wn / 2 + k][j - wn / 2 + l])*(x[i][j] - x[i - wn / 2 + k][j - wn / 2 + l]) / (lmdg*diat));
							N++;
						}
					}
				}
			}

			N = 0;
			s1 = 0; s2 = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						s1 += S[N] * x[i - wn / 2 + k][j - wn / 2 + l];
						s2 += S[N];
						N++;
					}
				}
			}
			y[i][j] = s1 / s2;
		}
	}
	delete[]S;
}


double similarity(double **x, int width, int height, int ci, int cj, int li, int lj, int s)
{
	int k, l;
	double sum, ik, pk;
	sum = 0;
	for (k = 0; k<s; k++)
	{
		for (l = 0; l<s; l++)
		{
			if (ci - s / 2 + k >= 0 && ci - s / 2 + k<width&&cj - s / 2 + l >= 0 && cj - s / 2 + l<height&&li - s / 2 + k >= 0 && li - s / 2 + k<width&&lj - s / 2 + l >= 0 && lj - s / 2 + l<height)
			{
				ik = x[ci - s / 2 + k][cj - s / 2 + l];
				pk = x[li - s / 2 + k][lj - s / 2 + l];
				if (ik == 0)
				{
					ik = 1;
				}
				if (pk == 0)
				{
					pk = 1;
				}
				sum += log(ik / pk + pk / ik);
			}
		}
	}
	return sum;
}


void non_local(double **x, double **&y, int width, int height, int r, int s, double h)
{
	int i, j, k, l, n;
	double *wip, zi;
	wip = new double[r*r];
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			n = 0; zi = 0;
			for (k = 0; k<r; k++)
			{
				for (l = 0; l<r; l++)
				{
					if (i - r / 2 + k >= 0 && i - r / 2 + k<width&&j - r / 2 + l >= 0 && j - r / 2 + l<height)
					{
						wip[n] = exp(-similarity(x, width, height, i, j, i - r / 2 + k, j - r / 2 + l, s) / h);
						zi += wip[n];
						n++;
					}
				}
			}
			n = 0;
			for (k = 0; k<r; k++)
			{
				for (l = 0; l<r; l++)
				{
					if (i - r / 2 + k >= 0 && i - r / 2 + k<width&&j - r / 2 + l >= 0 && j - r / 2 + l<height)
					{
						wip[n] = wip[n] / zi;
						n++;
					}
				}
			}

			n = 0; y[i][j] = 0;
			for (k = 0; k<r; k++)
			{
				for (l = 0; l<r; l++)
				{
					if (i - r / 2 + k >= 0 && i - r / 2 + k<width&&j - r / 2 + l >= 0 && j - r / 2 + l<height)
					{
						y[i][j] += wip[n] * x[i - r / 2 + k][j - r / 2 + l];
						n++;
					}
				}
			}
		}
	}
	delete[]wip;
}


double mlargest(double *x, int n, int m)
{
	int i, j, index;
	double s, a;
	if (m >= n)
	{
		s = 0;
		for (i = 0; i<n; i++)
		{
			s += x[i];
		}
		s = s / n;
	}
	else
	{
		s = 0;
		for (i = 0; i<m; i++)
		{
			a = x[i]; index = i;
			for (j = i + 1; j<n; j++)
			{
				if (x[j]>a)
				{
					a = x[j];
					index = j;
				}
			}
			x[index] = x[i];
			x[i] = a;
			s += x[i];
		}
		s = s / m;
	}
	return s;
}


void weighted_image(double **x, double **&y, int width, int height, int r)
{
	int i, j, k, l;
	double **local, **nlocal;

	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	double h = 10;
	local1(x, local, width, height, 3, 10, 10);
	non_local(x, nlocal, width, height, 21, 5, 2);

	double *ut, lmd;
	ut = new double[r*r];

	int n;

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			n = 0;
			for (k = 0; k<r; k++)
			{
				for (l = 0; l<r; l++)
				{
					if (i - r / 2 + k >= 0 && i - r / 2 + k<width&&j - r / 2 + l >= 0 && j - r / 2 + l<height)
					{
						ut[n] = exp(-similarity(x, width, height, i, j, i - r / 2 + k, j - r / 2 + l, r) / h);
						n++;
					}
				}
			}

			lmd = mlargest(ut, n, r - 1);
			//printf("%f\n",lmd);
			y[i][j] = (1 - lmd)*local[i][j] + lmd*nlocal[i][j];
		}
	}

	for (i = 0; i<width; i++)
	{
		delete[]local[i];
		delete[]nlocal[i];
	}
	delete[]local;
	delete[]nlocal;
	delete[]ut;
}


void cluster(double **x, double **&y, int width, int height, double m)
{
	int k, i, j;
	double v1, v2, v1new, v2new, s1, s2, s11, s12, s21, s22;
	double **u1, **u2;
	u1 = new double*[width];
	u2 = new double*[width];
	for (i = 0; i<width; i++)
	{
		u1[i] = new double[height];
		u2[i] = new double[height];
	}

	v1 = 255; v2 = 0; v1new = 0; v2new = 0;
	while ((v1 - v1new)*(v1 - v1new) + (v2 - v2new)*(v2 - v2new)>0.00001)
	{
		v1new = v1; v2new = v2;
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				if (x[i][j] == v1)
				{
					u1[i][j] = 1;
					u2[i][j] = 0;
				}
				else if (x[i][j] == v2)
				{
					u1[i][j] = 0;
					u2[i][j] = 1;
				}
				else
				{
					s1 = pow(abs(x[i][j] - v1), -2 / (m - 1));
					s2 = pow(abs(x[i][j] - v2), -2 / (m - 1));
					u1[i][j] = s1 / (s1 + s2);
					u2[i][j] = s2 / (s1 + s2);
				}
			}
		}
		s11 = 0; s12 = 0; s21 = 0; s22 = 0;
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				s11 += pow(u1[i][j], m)*x[i][j];
				s21 += pow(u2[i][j], m)*x[i][j];
				s12 += pow(u1[i][j], m);
				s22 += pow(u2[i][j], m);
			}
		}
		v1 = s11 / s12;
		v2 = s21 / s22;
		//printf("%f  %f\n",v1,v2);
	}

	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
		for (j = 0; j<height; j++)
		{
			if (u1[i][j]>u2[i][j])
			{
				y[i][j] = 0;
			}
			else
			{
				y[i][j] = 255;
			}
		}
	}
	for (i = 0; i<width; i++)
	{
		delete[]u1[i];
		delete[]u2[i];
	}
	delete[]u1;
	delete[]u2;
}

void fgfcm(double **x, double **&y, int width, int height)
{
	double **local;
	local1(x, local, width, height, 3, 10, 10);
	cluster(local, y, width, height, 2.0);
	int i;
	for (i = 0; i<width; i++)
	{
		delete[]local[i];
	}
	delete[]local;
}


void nlfcm(double**x, double **&y, int width, int height)
{
	double **nl;
	weighted_image(x, nl, width, height, 3);
	//cluster(nl, y, width, height, 2.0);
	//y = flicm(nl, width, height, 2.0, 3);
	y = nl;
	int i;
	for (i = 0; i<width; i++)
	{
		//delete[]nl[i];
	}
	//delete[]nl;
}





double min(double x1, double x2)
{
	if (x1>x2)
	{
		return x2;
	}
	else
	{
		return x1;
	}
}
void calgk(double **x, double **cu, double **u1, double **u2, double **g1, double **g2, double v1, double v2, int width, int height, int wn, int m)
{
	int i, j, k, l;
	double avrcu, n;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			avrcu = 0; n = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						avrcu += cu[i - wn / 2 + k][j - wn / 2 + l];
						n++;
					}
				}
			}
			avrcu = avrcu / n;
			if (cu[i][j]>avrcu)
			{
				g1[i][j] = 0; g2[i][j] = 0;
				for (k = 0; k<wn; k++)
				{
					for (l = 0; l<wn; l++)
					{
						if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
						{
							g1[i][j] += pow(1 - u1[i - wn / 2 + k][j - wn / 2 + l], m)*(x[i - wn / 2 + k][j - wn / 2 + l] - v1)*(x[i - wn / 2 + k][j - wn / 2 + l] - v1) / (2 + min((cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu)*(cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu), (avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])*(avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])));
							g2[i][j] += pow(1 - u2[i - wn / 2 + k][j - wn / 2 + l], m)*(x[i - wn / 2 + k][j - wn / 2 + l] - v2)*(x[i - wn / 2 + k][j - wn / 2 + l] - v2) / (2 + min((cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu)*(cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu), (avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])*(avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])));
						}
					}
				}
			}
			else
			{
				g1[i][j] = 0; g2[i][j] = 0;
				for (k = 0; k<wn; k++)
				{
					for (l = 0; l<wn; l++)
					{
						if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
						{
							g1[i][j] += pow(1 - u1[i - wn / 2 + k][j - wn / 2 + l], m)*(x[i - wn / 2 + k][j - wn / 2 + l] - v1)*(x[i - wn / 2 + k][j - wn / 2 + l] - v1) / (2 - min((cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu)*(cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu), (avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])*(avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])));
							g2[i][j] += pow(1 - u2[i - wn / 2 + k][j - wn / 2 + l], m)*(x[i - wn / 2 + k][j - wn / 2 + l] - v2)*(x[i - wn / 2 + k][j - wn / 2 + l] - v2) / (2 - min((cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu)*(cu[i - wn / 2 + k][j - wn / 2 + l] / avrcu), (avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])*(avrcu / cu[i - wn / 2 + k][j - wn / 2 + l])));
						}
					}
				}
			}
		}
	}
}
void rflicm(double **x, double **&y, int width, int height, int wn, int m)
{
	int i, j, k, l;
	double **cu;
	cu = new double*[width];
	for (i = 0; i<width; i++)
	{
		cu[i] = new double[height];
	}

	double u, diat, n;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			u = 0; n = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						u += x[i - wn / 2 + k][j - wn / 2 + l];
						n++;
					}
				}
			}
			u = u / n;
			diat = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						diat += (x[i - wn / 2 + k][j - wn / 2 + l] - u)*(x[i - wn / 2 + k][j - wn / 2 + l]);
					}
				}
			}
			diat = diat / n;
			if (u == 0)
			{
				cu[i][j] = 0;
			}
			else
			{
				cu[i][j] = diat / (u*u);
			}
		}
	}
	double **g1, **g2, **u1, **u2;
	g1 = new double*[width];
	g2 = new double*[width];
	u1 = new double*[width];
	u2 = new double*[width];
	for (i = 0; i<width; i++)
	{
		g1[i] = new double[height];
		g2[i] = new double[height];
		u1[i] = new double[height];
		u2[i] = new double[height];
		for (j = 0; j<height; j++)
		{
			u1[i][j] = 0; u2[i][j] = 0;
		}
	}
	double v1, v2, v1new, v2new, mid1, mid2, s11, s12, s21, s22;
	v1 = 0; v2 = 255; v1new = 0; v2new = 0;
	while ((v1 - v1new)*(v1 - v1new) + (v2 - v2new)*(v2 - v2new)>0.0001)
	{
		v1new = v1; v2new = v2;
		calgk(x, cu, u1, u2, g1, g2, v1, v2, width, height, wn, m);

		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				mid1 = (x[i][j] - v1)*(x[i][j] - v1) + g1[i][j];
				mid2 = (x[i][j] - v2)*(x[i][j] - v2) + g2[i][j];
				u1[i][j] = 1 / (1 + pow(mid1 / mid2, 1 / (double)(m - 1)));
				u2[i][j] = 1 / (1 + pow(mid2 / mid1, 1 / (double)(m - 1)));
			}
		}

		s11 = 0; s12 = 0; s21 = 0; s22 = 0;
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				s11 += pow(u1[i][j], m)*x[i][j];
				s12 += pow(u1[i][j], m);
				s21 += pow(u2[i][j], m)*x[i][j];
				s22 += pow(u2[i][j], m);
			}
		}
		v1 = s11 / s12;
		v2 = s21 / s22;
	}

	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (u1[i][j]>u2[i][j])
			{
				y[i][j] = 0;
			}
			else
			{
				y[i][j] = 255;
			}
		}
		delete[]cu[i];
		delete[]g1[i];
		delete[]g2[i];
		delete[]u1[i];
		delete[]u2[i];
	}
	delete[]cu;
	delete[]g1;
	delete[]g2;
	delete[]u1;
	delete[]u2;
}
