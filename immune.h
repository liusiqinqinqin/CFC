/******************************************************************************************************************************************************
code for artificial immune system method which is used to classify the difference image.

1main function:
immune: code for classifying the difference image by artificial immune system.

input:
x1: a matrix, the first image with the data type of double;
x2: a matrix, the second image with the data type of double;
width: the width of the two image;
height: the height of the two image;
num_ab: the number of antibodies;
num_mc: the number of memory cells;
num_select: the number of selected antobodies during each iteration;
clonal_rate: the number of multiple cloned antibodies;
num_replace: the number of replaced antibodies;
wn: the size of neighborhood (wn*wn).
num_ab, num_mc, clonal_rate, num_replace are the parameters used in the artificial immune system, the reference values of them are:
num_ab=20;
num_mc=3;
num_select=5;
clonal_rate=20;
num_replace=5;

output:
the return of the function: a matrix, the changed map with the data type of double.

The rest functions are intermediate functions which are steps of artificial immune system.
******************************************************************************************************************************************************/
#include"flicm.h"

#define e 2.718281828459

double f_KI(double **x, long width, long height, int sd)
{
	int i, j;
	double u1, u2, p1, p2, diat1, diat2, n1, n2;                 //1 change,    2 unchange
	u1 = 0; u2 = 0; n1 = 0; n2 = 0;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (x[i][j] <= (double)sd)
			{
				u1 += x[i][j];
				n1++;
			}
			else
			{
				u2 += x[i][j];
				n2++;
			}
		}
	}
	if (n1 != 0)
	{
		u1 = u1 / n1;
	}
	if (n2 != 0)
	{
		u2 = u2 / n2;
	}
	p1 = n1 / (n1 + n2);
	p2 = n2 / (n1 + n2);
	diat1 = 0; diat2 = 0;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (x[i][j] <= (double)sd)
			{
				diat1 += (x[i][j] - u1)*(x[i][j] - u1);
			}
			else
			{
				diat2 += (x[i][j] - u2)*(x[i][j] - u2);
			}
		}
	}
	if (n1 != 0)
	{
		diat1 = sqrt(diat1 / n1);
	}
	if (n2 != 0)
	{
		diat2 = sqrt(diat2 / n2);
	}
	if (diat1 == 0 || diat2 == 0 || p1 == 0 || p2 == 0)
	{
		return 0;
	}
	else
	{
		return 1 + 2 * (p1*log(diat1) / log(e) + p2*log(diat2) / log(e)) - 4 * (p1*log(p1) / log(e) + p2*log(p2) / log(e));
		//return -p1*log(p1)-p2*log(p2);
	}
	//return 1+2*(p1*log(diat1+1)/log(e)+p2*log(diat2+1)/log(e))-4*(p1*log(p1+1)/log(e)+p2*log(p2+1)/log(e));

}
int threshold(double **x, long width, long height)
{
	int i, index;
	double min, fun;
	min = f_KI(x, width, height, 0);
	index = 0;
	for (i = 1; i<256; i++)
	{
		fun = f_KI(x, width, height, i);
		if (min<fun)
		{
			min = fun;
			index = i;
		}
	}
	return index + 30;
}


double distance(double x1, double x2)
{
	if (x1>x2)
	{
		return x1 - x2;
	}
	else
	{
		return x2 - x1;
	}
}

double ab_ag_distance(double ab, double **x, int width, int height, int pos_i, int pos_j, double **wij, int wn)
{
	double s, num;
	int i, j, px, py;
	s = 0; num = 0;
	for (i = 0; i<wn; i++)
	{
		for (j = 0; j<wn; j++)
		{
			px = pos_i - wn / 2 + i;
			py = pos_j - wn / 2 + j;
			if (px >= 0 && px<width&&py >= 0 && py<height)
			{
				s += wij[i][j] * (x[px][py] - ab)*(x[px][py] - ab);
				num += wij[i][j];
			}
		}
	}

	s = s / num;
	return s;
}

void select(double *x, double *x_distance, int num_x, double *y, double *y_distance, int num_y)
{
	int i, j;
	int *label;
	label = new int[num_x];
	for (i = 0; i<num_x; i++)
	{
		label[i] = 0;
	}

	double min;
	int index;
	for (i = 0; i<num_y; i++)
	{
		min = 65530.0;
		for (j = 0; j<num_x; j++)
		{
			if (x_distance[j]<min&&label[j] == 0)
			{
				min = x_distance[j];
				index = j;
			}
		}
		y[i] = x[index];
		y_distance[i] = x_distance[index];
		label[index] = 1;
	}

	delete[]label;
}


double Nrandom(void)
{
	double x;
	int i;
	x = 0;
	for (i = 0; i<12; i++)
	{
		x += (double)rand() / (double)RAND_MAX;
	}
	return x - 6;
}
void mutate(double *selected, double *select_distance, int num_select, int clonal_rate, double *replace, int num_replace, double **x, int width, int height, int pos_i, int pos_j, double **wij, int wn)
{
	double *MU, *mu_distance;
	int i, j, num_mu;
	num_mu = num_select*clonal_rate;
	MU = new double[num_mu];
	mu_distance = new double[num_mu];

	int index;
	index = 0;
	for (i = 0; i<num_select; i++)
	{
		for (j = 0; j<clonal_rate; j++)
		{
			MU[index] = selected[i] + Nrandom();
			index++;
		}
	}

	for (i = 0; i<num_mu; i++)
	{
		mu_distance[i] = ab_ag_distance(MU[i], x, width, height, pos_i, pos_j, wij, wn);
	}

	double *replace_distance;
	replace_distance = new double[num_replace];
	select(MU, mu_distance, num_mu, replace, replace_distance, num_replace);

	delete[]MU;
	delete[]mu_distance;
	delete[]replace_distance;
}

void replace_AB(double *AB, double *ab_distance, int num_ab, double *replace, int num_replace)
{
	int i, j;
	int *label;
	label = new int[num_ab];
	for (i = 0; i<num_ab; i++)
	{
		label[i] = 0;
	}

	double max;
	int index;
	for (i = 0; i<num_replace; i++)
	{
		max = 0;
		for (j = 0; j<num_ab; j++)
		{
			if (ab_distance[j]>max&&label[j] == 0)
			{
				max = ab_distance[j];
				index = j;
			}
		}
		AB[index] = replace[i];
		label[index] = 1;
	}
	delete[]label;
}



double update_AB(double *AB, int num_ab, double **x, int width, int height, int pos_i, int pos_j, double **wij, int wn, int clonal_rate, int num_select, int num_replace)
{
	int i, j;
	double *ab_distance;
	ab_distance = new double[num_ab];

	for (i = 0; i<num_ab; i++)
	{
		ab_distance[i] = ab_ag_distance(AB[i], x, width, height, pos_i, pos_j, wij, wn);
	}

	double *selected, *select_distance;
	selected = new double[num_select];
	select_distance = new double[num_select];
	select(AB, ab_distance, num_ab, selected, select_distance, num_select);

	double *replace;
	replace = new double[num_replace];
	mutate(selected, select_distance, num_select, clonal_rate, replace, num_replace, x, width, height, pos_i, pos_j, wij, wn);

	replace_AB(AB, ab_distance, num_ab, replace, num_replace);
	double candidate;
	candidate = replace[0];

	delete[]ab_distance;
	delete[]selected;
	delete[]select_distance;
	delete[]replace;

	return candidate;
}

void update_MC(double *MC, int n_mc, double candidate, double **x, int width, int height, int pos_i, int pos_j, double**wij, int wn)
{
	int i;
	double min, dis;
	int index;
	min = distance(candidate, MC[0]);
	index = 0;
	for (i = 1; i<n_mc; i++)
	{
		dis = distance(candidate, MC[i]);
		if (dis<min)
		{
			min = dis;
			index = i;
		}
	}

	if (ab_ag_distance(candidate, x, width, height, pos_i, pos_j, wij, wn)<ab_ag_distance(MC[index], x, width, height, pos_i, pos_j, wij, wn))
	{
		MC[index] = candidate;
	}
}


void recognize(double *MC0, int num_mc0, double *MC1, int num_mc1, double **x, int width, int height, int pos_i, int pos_j, double **wij, int wn, int &clas, double &blong0, double &blong1)
{
	int i;
	double min0, min1, distance;
	min0 = ab_ag_distance(MC0[0], x, width, height, pos_i, pos_j, wij, wn);
	for (i = 1; i<num_mc0; i++)
	{
		distance = ab_ag_distance(MC0[i], x, width, height, pos_i, pos_j, wij, wn);
		if (distance<min0)
		{
			min0 = distance;
		}
	}

	min1 = ab_ag_distance(MC1[0], x, width, height, pos_i, pos_j, wij, wn);
	for (i = 1; i<num_mc1; i++)
	{
		distance = ab_ag_distance(MC1[i], x, width, height, pos_i, pos_j, wij, wn);
		if (distance<min1)
		{
			min1 = distance;
		}
	}

	if (min0<min1)
	{
		clas = 0;
	}
	else
	{
		clas = 1;
	}

	blong0 = min1 / (min0 + min1);
	blong1 = min0 / (min0 + min1);
}

void initial_AB(double *&AB, int num_ab, int clas, int T, double **x, int width, int height)
{
	int i, j, index, px, py;
	AB = new double[num_ab];


	for (i = 0; i<num_ab; i++)
	{
		if (clas == 0)
		{
			index = 1;
			while (index == 1)
			{
				px = rand() % width;
				py = rand() % height;
				index = 0;
				for (j = 0; j<i; j++)
				{
					if (distance(x[px][py], AB[j]) == 0)
					{
						index = 1;
					}
				}
				if (x[px][py]>T)
				{
					index = 1;
				}
			}

			AB[i] = x[px][py];
		}
		else
		{
			index = 1;
			while (index == 1)
			{
				px = rand() % width;
				py = rand() % height;
				index = 0;
				for (j = 0; j<i; j++)
				{
					if (distance(x[px][py], AB[j]) == 0)
					{
						index = 1;
					}
				}
				if (x[px][py]<T)
				{
					index = 1;
				}
			}

			AB[i] = x[px][py];
		}
	}
}
void initial_MC(double *AB, int num_ab, double *&MC, int num_mc)
{
	double mean;
	int i, j;

	mean = 0;

	for (i = 0; i<num_ab; i++)
	{
		mean += AB[i];
	}
	mean = mean / (double)num_ab;


	double min;
	int index;
	min = distance(AB[0], mean);
	index = 0;
	for (i = 0; i<num_mc; i++)
	{
		if (distance(AB[i], mean)<min)
		{
			min = distance(AB[i], mean);
			index = i;
		}
	}

	MC = new double[num_mc];


	MC[0] = AB[index];


	int *label;
	label = new int[num_ab];
	for (i = 0; i<num_ab; i++)
	{
		label[i] = 0;
	}
	label[index] = 1;

	int k, index1, l;
	double max, s;
	for (k = 1; k<num_mc; k++)
	{
		max = -1;
		for (i = 0; i<num_ab; i++)
		{
			if (label[i] == 0)
			{
				s = 0;
				for (j = 0; j<num_ab; j++)
				{
					if (j != i)
					{
						min = distance(MC[0], AB[j]);
						index = 0;
						for (l = 1; l<k; l++)
						{
							if (distance(MC[l], AB[j])<min)
							{
								min = distance(MC[l], AB[j]);
								index = l;
							}
						}
						if (min - distance(AB[i], AB[j])>0)
						{
							s += min - distance(AB[i], AB[j]);
						}
					}
				}
				if (s>max)
				{
					max = s;
					index1 = i;
				}
			}
		}
		MC[k] = AB[index1];
		label[index1] = 1;
	}
	delete[]label;
}

/*void cal_weight(double ****&x_weight, double **x, int width, int height, int wn)
{
	int i, j, k, l;
	x_weight = new double***[width];
	for (i = 0; i<width; i++)
	{
		x_weight[i] = new double**[height];
		for (j = 0; j<height; j++)
		{
			x_weight[i][j] = new double*[wn];
			for (k = 0; k<wn; k++)
			{
				x_weight[i][j][k] = new double[wn];
			}
		}
	}

	int px, py;
	double **C;
	C = new double*[width];
	for (i = 0; i<width; i++)
	{
		C[i] = new double[height];
	}

	double var, avr, num;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			avr = 0; var = 0; num = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						avr += x[px][py];
						num++;
					}
				}
			}
			avr = avr / num;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						var += (x[px][py] - avr)*(x[px][py] - avr);
					}
				}
			}
			var = var / num;
			C[i][j] = var / (avr*avr);
		}
	}

	double **avrC;
	avrC = new double*[width];
	for (i = 0; i<width; i++)
	{
		avrC[i] = new double[height];
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			avrC[i][j] = 0;
			num = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						avrC[i][j] += C[px][py];
						num += 1;
					}
				}
			}
			avrC[i][j] = avrC[i][j] / num;
		}
	}

	double **keci, **yita, sum;
	keci = new double*[wn];
	yita = new double*[wn];
	for (i = 0; i<wn; i++)
	{
		keci[i] = new double[wn];
		yita[i] = new double[wn];
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			sum = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						keci[k][l] = exp(-(C[px][py] - avrC[i][j]));
						sum += keci[k][l];
					}
				}
			}
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						yita[k][l] = keci[k][l] / sum;
					}
				}
			}
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					px = i - wn / 2 + k;
					py = j - wn / 2 + l;
					if (px >= 0 && px<width&&py >= 0 && py<height)
					{
						if (C[px][py]<avrC[i][j])
						{
							x_weight[i][j][k][l] = (2 + yita[k][l]) / (sqrt((double)((i - px)*(i - px) + (j - py)*(j - py))) + 1);
						}
						else
						{
							x_weight[i][j][k][l] = (2 - yita[k][l]) / (sqrt((double)((i - px)*(i - px) + (j - py)*(j - py))) + 1);
						}
					}
				}
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		delete[]C[i];
		delete[]avrC[i];
	}
	for (i = 0; i<wn; i++)
	{
		delete[]keci[i];
		delete[]yita[i];
	}
	delete[]C;
	delete[]avrC;
	delete[]keci;
	delete[]yita;

}*/
void cal_wij(double **wij, int wn, double ****x_weight, double **x, double ***map, int width, int height, int pos_i, int pos_j)
{
	int i, j, px, py;
	for (i = 0; i<wn; i++)
	{
		for (j = 0; j<wn; j++)
		{
			px = pos_i - wn / 2 + i;
			py = pos_j - wn / 2 + j;
			if (px >= 0 && px<width&&py >= 0 && py<height)//&&(px!=pos_i||py!=pos_j))
			{
				wij[i][j] = x_weight[pos_i][pos_j][i][j];

				/*if(map[px][py][0]!=map[pos_i][pos_j][0])
				{
				wij[i][j]=x_weight[pos_i][pos_j][i][j]*(1-map[px][py][(int)map[pos_i][pos_j][0]+1]);
				}
				else
				{
				wij[i][j]=x_weight[pos_i][pos_j][i][j]*(1-map[px][py][(int)map[pos_i][pos_j][0]+1]);
				}*/

			}
		}
	}
	wij[(int)(wn / 2)][(int)(wn / 2)] = 1;
}


double** immune(double **x, int width, int height, int num_ab, int num_mc, int num_select, int clonal_rate, int num_replace, int wn)
{
	int i, j, k, l;
	double *MC0, *MC1, *AB0, *AB1;
	int T, num_mc0, num_mc1;
	num_mc0 = num_mc; num_mc1 = num_mc;
	//num_mc0=7;num_mc1=3;
	T = threshold(x, width, height);
	initial_AB(AB0, num_ab, 0, T, x, width, height);
	initial_AB(AB1, num_ab, 1, T, x, width, height);
	initial_MC(AB0, num_ab, MC0, num_mc0);
	initial_MC(AB1, num_ab, MC1, num_mc1);

	double ****x_weight, **wij;
	wij = new double*[wn];
	for (i = 0; i<wn; i++)
	{
		wij[i] = new double[wn];
	}

	double ***map, ***map1;
	map = new double**[width];
	map1 = new double**[width];
	for (i = 0; i<width; i++)
	{
		map[i] = new double*[height];
		map1[i] = new double*[height];
		for (j = 0; j<height; j++)
		{
			map[i][j] = new double[3];
			map1[i][j] = new double[3];
			if (x[i][j]>(double)T)
			{
				map[i][j][0] = 1;
				map[i][j][1] = 0;
				map[i][j][2] = 1;
			}
			else
			{
				map[i][j][0] = 0;
				map[i][j][1] = 1;
				map[i][j][2] = 0;
			}
			for (k = 0; k<3; k++)
			{
				map1[i][j][k] = map[i][j][k];
			}

		}
	}

	cal_weight(x_weight, x, width, height, wn);


	int clas;
	double u0, u1, candidate;
	for (k = 0; k<10; k++)
	{
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				cal_wij(wij, wn, x_weight, x, map1, width, height, i, j);
				recognize(MC0, num_mc0, MC1, num_mc1, x, width, height, i, j, wij, wn, clas, u0, u1);
				if (clas == 0)
				{
					candidate = update_AB(AB0, num_ab, x, width, height, i, j, wij, wn, clonal_rate, num_select, num_replace);
					update_MC(MC0, num_mc, candidate, x, width, height, i, j, wij, wn);
				}
				else
				{
					candidate = update_AB(AB1, num_ab, x, width, height, i, j, wij, wn, clonal_rate, num_select, num_replace);
					update_MC(MC1, num_mc, candidate, x, width, height, i, j, wij, wn);
				}
				map[i][j][0] = (double)clas;
				map[i][j][1] = u0;
				map[i][j][2] = u1;
			}

		}
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				for (l = 0; l<3; l++)
				{
					map1[i][j][l] = map[i][j][l];
				}
			}
		}

		printf("%d\n", k);
	}

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
			cal_wij(wij, wn, x_weight, x, map1, width, height, i, j);
			recognize(MC0, num_mc0, MC1, num_mc1, x, width, height, i, j, wij, wn, clas, u0, u1);
			y[i][j] = (double)clas*255.0;
			/*map1[i][j][0]=(double)clas;
			map1[i][j][1]=u0;
			map1[i][j][2]=u1;*/
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			for (k = 0; k<wn; k++)
			{
				delete[]x_weight[i][j][k];
			}
			delete[]x_weight[i][j];
			delete[]map[i][j];
		}
		delete[]x_weight[i];
		delete[]map[i];
	}
	for (i = 0; i<wn; i++)
	{
		delete[]wij[i];
	}
	delete[]wij;
	delete[]AB0;
	delete[]MC0;
	delete[]AB1;
	delete[]MC1;

	return y;
}