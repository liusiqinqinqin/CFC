
#include"flicm.h"



double** imageread(char *s, int &width, int &height)
{
	int i, j;
	double **y;
	IplImage *image;

	image = cvLoadImage(s);

	width = image->width;
	height = image->height;

	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	CvScalar x;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			x = cvGet2D(image, j, i);
			y[i][j] = (double)x.val[0];
		}
	}

	cvReleaseImage(&image);
	return y;
}
void imagewrite(double **x, int width, int height, char*s, char *name)
{
	int i, j;
	IplImage *image;
	image = cvLoadImage(s);
	CvScalar mid;

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			mid.val[0] = (int)(x[i][j] + 0.5);
			mid.val[1] = (int)(x[i][j] + 0.5);
			mid.val[2] = (int)(x[i][j] + 0.5);
			cvSet2D(image, j, i, mid);
		}
	}

	char save[100];
	i = 0;
	while (s[i] != '\0')
	{
		i++;
	}
	while (s[i] != 92)
	{
		i--;
	}
	for (j = 0; j <= i; j++)
	{
		save[j] = s[j];
	}
	for (j = 0; j<strlen(name) + 1; j++)
	{
		save[i + 1 + j] = name[j];
	}

	cvSaveImage(save, image);

	cvReleaseImage(&image);
}




void imageconvert(double **x, int width, int height)
{
	int i, j;
	double max, min;
	max = -1000; min = 1000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (x[i][j]>max)
			{
				max = x[i][j];
			}
			if (x[i][j]<min)
			{
				min = x[i][j];
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			x[i][j] = 255 * (x[i][j] - min) / (max - min);
		}
	}
}

/*void imageread(char *s,double **&x,int &width,int &height)
{

IplImage *image;
image=cvLoadImage(s);

int i,j;
width=image->width;
height=image->height;
x=new double*[width];
for(i=0;i<width;i++)
{
x[i]=new double[height];
}

int mid;
if(image->nChannels==1)
{
for(i=0;i<height;i++)
{
for(j=0;j<width;j++)
{
image->imageData[i*width+j]+=128;
mid=image->imageData[i*width+j]+128;
image->imageData[i*width+j]-=128;
x[j][i]=(double)mid;
}
}
}
if(image->nChannels==3)
{
for(i=0;i<height;i++)
{
for(j=0;j<width;j++)
{
image->imageData[(i*width+j)*3]+=128;
mid=image->imageData[(i*width+j)*3]+128;
image->imageData[(i*width+j)*3]-=128;
x[j][i]=(double)mid;
}
}
}

cvReleaseImage(&image);
}*/

void NR(double **x1, double **x2, int width, int height, double **&y, int wn)
{
	int i, j, k, l;
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}
	double num, u, diat, min, max, smin, smax, sp, np, emax, emin;
	emax = -1000; emin = 1000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			u = 0;
			num = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						u += x1[i - wn / 2 + k][j - wn / 2 + l];
						u += x2[i - wn / 2 + k][j - wn / 2 + l];
						num += 2;
					}
				}
			}
			u = u / num;

			diat = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						diat += (x1[i - wn / 2 + k][j - wn / 2 + l] - u)*(x1[i - wn / 2 + k][j - wn / 2 + l] - u);
						diat += (x2[i - wn / 2 + k][j - wn / 2 + l] - u)*(x2[i - wn / 2 + k][j - wn / 2 + l] - u);
					}
				}
			}
			diat = sqrt(diat / num);
			//diat/=10;
			diat = diat / u;

			smin = 0; smax = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						if (i - wn / 2 + k != i || j - wn / 2 + l)
						{
							if (x1[i - wn / 2 + k][j - wn / 2 + l]>x2[i - wn / 2 + k][j - wn / 2 + l])
							{
								smin += x2[i - wn / 2 + k][j - wn / 2 + l];
								smax += x1[i - wn / 2 + k][j - wn / 2 + l];
							}
							else
							{
								smin += x1[i - wn / 2 + k][j - wn / 2 + l];
								smax += x2[i - wn / 2 + k][j - wn / 2 + l];
							}
						}
					}
				}
			}

			if (x1[i][j]>x2[i][j])
			{
				min = x2[i][j];
				max = x1[i][j];
			}
			else
			{
				min = x1[i][j];
				max = x2[i][j];
			}
			if (max == 0)
			{
				sp = 1;
			}
			else
			{
				sp = min / max;
			}

			if (smax == 0)
			{
				np = 1;
			}
			else
			{
				np = smin / smax;
			}

			y[i][j] = diat*sp + (1 - diat)*np;
			if (y[i][j]>emax)
			{
				emax = y[i][j];
			}
			if (y[i][j]<emin)
			{
				emin = y[i][j];
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 255 - 255 * (y[i][j] - emin) / (emax - emin);
		}
	}
}
void sub(double **x1, double **x2, int width, int height, double **&y)
{
	int i, j;
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	double min, max;
	min = 1000; max = -1000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = abs(x1[i][j] - x2[i][j]);
			if (y[i][j]>max)
			{
				max = y[i][j];
			}
			if (y[i][j]<min)
			{
				min = y[i][j];
			}
		}
	}
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 255 * (y[i][j] - min) / (max - min);
		}
	}
}

void log_ratio(double **x1, double **x2, int width, int height, double **&y)
{
	int i, j;
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	double min, max;
	min = 1000; max = -1000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			/*if(x1[i][j]==0||x2[i][j]==0)
			{
			y[i][j]=255;
			}
			else
			{
			y[i][j]=abs(log(x1[i][j])-log(x2[i][j]));
			}*/
			y[i][j] = log(x1[i][j] + 1) - log(x2[i][j] + 1);
			if (y[i][j]>max)
			{
				max = y[i][j];
			}
			if (y[i][j]<min)
			{
				min = y[i][j];
			}
		}
	}

	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 255.0*(y[i][j] - min) / (max - min);
		}
	}
}

void mean_ratio(double **x1, double **x2, int width, int height, double **&y, int wn)
{
	int i, j, k, l;
	y = new double*[width];
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	double u1, u2, max, min, num;
	max = -1000; min = 1000;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			u1 = 0; u2 = 0; num = 0;
			for (k = 0; k<wn; k++)
			{
				for (l = 0; l<wn; l++)
				{
					if (i - wn / 2 + k >= 0 && i - wn / 2 + k<width&&j - wn / 2 + l >= 0 && j - wn / 2 + l<height)
					{
						u1 += x1[i - wn / 2 + k][j - wn / 2 + l];
						u2 += x2[i - wn / 2 + k][j - wn / 2 + l];
						num++;
					}
				}
			}
			u1 /= num;
			u2 /= num;
			/*if(u1>u2)
			{
			if(u1!=0)
			{
			y[i][j]=u2/u1;
			}
			else
			{
			y[i][j]=1;
			}
			}
			else
			{
			if(u2!=0)
			{
			y[i][j]=u1/u2;
			}
			else
			{
			y[i][j]=1;
			}
			}*/
			y[i][j] = log(u1 + 0.0001) - log(u2 + 0.0001);


			if (y[i][j]>max)
			{
				max = y[i][j];
			}
			if (y[i][j]<min)
			{
				min = y[i][j];
			}
		}
	}


	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			y[i][j] = 255 - 255 * (y[i][j] - min) / (max - min);
		}
	}
}

void nsct1(double **x1, double **x2, int width, int height, double **&y)
{
	y = new double*[width];
	int i;
	for (i = 0; i<width; i++)
	{
		y[i] = new double[height];
	}

	double **lr, **mr;
	log_ratio(x1, x2, width, height, lr);
	mean_ratio(x1, x2, width, height, mr, 3); printf("ok\n");
	nsct_fusion(lr, mr, y, width, height);

	imageconvert(y, width, height);

	for (i = 0; i<width; i++)
	{
		delete[]lr[i];
		delete[]mr[i];
	}
	delete[]lr;
	delete[]mr;
}







void roc(double **x, double **gt, long width, long height, char *fname)
{
	int i, j, k;
	int re;

	int *TP, *TN, *FP, *FN;
	TP = (int*)malloc(256 * sizeof(int));
	TN = (int*)malloc(256 * sizeof(int));
	FP = (int*)malloc(256 * sizeof(int));
	FN = (int*)malloc(256 * sizeof(int));
	for (k = 0; k<256; k++)
	{
		TP[k] = 0; TN[k] = 0; FP[k] = 0; FN[k] = 0;
		for (i = 0; i<width; i++)
		{
			for (j = 0; j<height; j++)
			{
				if (x[i][j]<k)
				{
					re = 0;
				}
				else
				{
					re = 255;
				}
				if (re == 0 && gt[i][j] == 0)
				{
					TP[k]++;
				}
				if (re == 255 && gt[i][j] == 255)
				{
					TN[k]++;
				}
				if (re == 0 && gt[i][j] == 255)
				{
					FP[k]++;
				}
				if (re == 255 && gt[i][j] == 0)
				{
					FN[k]++;
				}
			}
		}
	}

	double *TPR, *FPR;
	TPR = (double*)malloc(256 * sizeof(double));
	FPR = (double*)malloc(256 * sizeof(double));
	for (i = 0; i<256; i++)
	{
		TPR[i] = (double)TP[i] / (double)(TP[i] + FN[i]);
		FPR[i] = (double)FP[i] / (double)(TN[i] + FP[i]);
	}



	double s = 0;
	for (i = 1; i<256; i++)
	{
		s += (TPR[i] + TPR[i - 1])*(FPR[i] - FPR[i - 1]) / 2;
	}
	printf("roc曲线面积：%f\n", s);


	FILE *plot;
	fopen_s(&plot,fname, "w+");
	for (i = 0; i<256; i++)
	{
		fprintf(plot, "%f ", TPR[i]);
	}
	for (i = 0; i<256; i++)
	{
		fprintf(plot, "%f ", FPR[i]);
	}
	fclose(plot);

	free(TN); free(FN); free(TP); free(FP); free(TPR); free(FPR);
}

void statistic(double **result, double**gt, int width, int height)
{
	int i, j;
	int TP, TN, FP, FN;
	TP = 0; TN = 0; FP = 0; FN = 0;
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (result[i][j] == 0 && gt[i][j] == 0)
			{
				TP++;
			}
			if (result[i][j] == 255 && gt[i][j] == 0)
			{
				TN++;
			}
			if (result[i][j] == 0 && gt[i][j] == 255)
			{
				FP++;
			}
			if (result[i][j] == 255 && gt[i][j] == 0)
			{
				FN++;
			}
		}
	}

	printf("FN=%d   FP=%d\noe=%d\nPCC=%f\n", FN, FP, FN + FP, (double)(TP + TN) / (double)(TP + TN + FP + FN));
}


void log(char *s1, char *s2, char *s3)
{
	double **x1, **x2, **y;
	int width, height, i, j;

	x1 = imageread(s1, width, height);
	printf("image1 ready\n");


	x2 = imageread(s2, width, height);
	printf("image2 ready\n");

	log_ratio(x1, x2, width, height, y);
	//imagewrite(y,width,height,s1,"log_ratio1.bmp");


	imagewrite(y, width, height, s1, "log.bmp");

	double**fg;
	//local1(y,fg,width,height,5,0.1,0.2);
	//non_local(y,fg,width,height,3,7,2);
	//weighted_image(y,fg,width,height,7);
	cluster(y, fg, width, height, 2);
	//fgfcm(y,fg,width,height);
	//nlfcm(y,fg,width,height);
	//rflicm(y,fg,width,height,5,2);


	double **gt;

	gt = imageread(s3, width, height);
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (gt[i][j]<128)
			{
				gt[i][j] = 0;
			}
			else
			{
				gt[i][j] = 255;
			}
		}
	}

	//imagedisplay(gt,width,height);

	roc(y, gt, width, height, "log_ratio.txt");
	statistic(fg, gt, width, height);

	printf("log_ratio complete\n\n");
	getchar();

	for (i = 0; i<width; i++)
	{
		delete[] x1[i];
		delete[] x2[i];
		delete[] y[i];
		delete[] gt[i];
		delete[]fg[i];
	}
	delete[]x1;
	delete[]x2;
	delete[]y;
	delete[]gt;
	delete[]fg;
}
void mean(char *s1, char *s2, char *s3)
{
	double **x1, **x2, **y;
	int width, height, i, j;

	x1 = imageread(s1, width, height);
	printf("image1 ready\n");

	x2 = imageread(s2, width, height);
	printf("image2 ready\n");

	mean_ratio(x1, x2, width, height, y, 7);



	imagewrite(y, width, height, s1, "mean-ratio.bmp");
	double**fg;
	//local1(y,fg,width,height,5,0.1,0.2);
	//non_local(y,fg,width,height,3,7,2);
	//weighted_image(y,fg,width,height,7);
	cluster(y, fg, width, height, 2);
	//fgfcm(y,fg,width,height);
	//nlfcm(y,fg,width,height);
	//rflicm(y,fg,width,height,5,2);


	double **gt;

	gt = imageread(s3, width, height);
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (gt[i][j]<128)
			{
				gt[i][j] = 0;
			}
			else
			{
				gt[i][j] = 255;
			}
		}
	}

	//imagedisplay(gt,width,height);

	roc(y, gt, width, height, "mean_ratio.txt");
	statistic(fg, gt, width, height);

	printf("mean_ratio complete\n\n");
	getchar();

	for (i = 0; i<width; i++)
	{
		delete[] x1[i];
		delete[] x2[i];
		delete[] y[i];
		delete[] gt[i];
		delete[]fg[i];
	}
	delete[]x1;
	delete[]x2;
	delete[]y;
	delete[]gt;
	delete[]fg;
}

void ct(char *s1, char *s2, char *s3)
{
	double **x1, **x2, **y;
	int width, height, i, j;

	x1 = imageread(s1, width, height);
	printf("image1 ready\n");


	x2 = imageread(s2, width, height);
	printf("image2 ready\n");


	nsct1(x1, x2, width, height, y);
	//log_ratio(x1, x2, width, height, y);


	imagewrite(y, width, height, s1, "contourlet_fusion.bmp");

	double**fg;
	//local1(y,fg,width,height,5,0.1,0.2);
	//non_local(y,fg,width,height,3,7,2);
	//weighted_image(y,fg,width,height,7);
	//cluster(y, fg, width, height, 2);
	//fgfcm(y,fg,width,height);
	nlfcm(y,fg,width,height);
	//rflicm(y,fg,width,height,5,2);
	imagewrite(fg, width, height, s1, "contourlet_fusion_result.bmp");

	double **gt;

	gt = imageread(s3, width, height);
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (gt[i][j]<128)
			{
				gt[i][j] = 0;
			}
			else
			{
				gt[i][j] = 255;
			}
		}
	}

	//imagedisplay(gt,width,height);

	roc(y, gt, width, height, "nsct.txt");
	statistic(fg, gt, width, height);

	printf("nsct fusion complete\n\n");
	//getchar();

	for (i = 0; i<width; i++)
	{
		delete[] x1[i];
		delete[] x2[i];
		delete[] y[i];
		delete[] gt[i];
		delete[]fg[i];
	}
	delete[]x1;
	delete[]x2;
	delete[]y;
	delete[]gt;
	delete[]fg;
}

void comparafcm(char *s1, char *s2, char *s3)
{
	double **x1, **x2, **y;
	int width, height, i, j;

	x1 = imageread(s1, width, height);
	printf("image1 ready\n");
	x2 = imageread(s2, width, height);
	printf("image2 ready\n");
	nsct1(x1, x2, width, height, y);
	//wavelet(x1,x2,width,height,y);
	//log_ratio(x1,x2,width,height,y);



	double **gt;
	gt = imageread(s3, width, height);
	for (i = 0; i<width; i++)
	{
		for (j = 0; j<height; j++)
		{
			if (gt[i][j]<128)
			{
				gt[i][j] = 0;
			}
			else
			{
				gt[i][j] = 255;
			}
		}
	}

	double **fg;
	printf("start fcm\n");
	cluster(y, fg, width, height, 2.0);

	statistic(fg, gt, width, height);
	for (i = 0; i<width; i++)
	{
		delete[]fg[i];
	}
	delete[]fg;
	printf("fcm complete\n\n");
	getchar();


	printf("start rflicm\n");
	rflicm(y, fg, width, height, 3, 2);

	statistic(fg, gt, width, height);
	for (i = 0; i<width; i++)
	{
		delete[]fg[i];
	}
	delete[]fg;
	printf("rflifcm complete\n\n");
	getchar();



	printf("start fgfcm\n");
	fgfcm(y, fg, width, height);

	statistic(fg, gt, width, height);
	for (i = 0; i<width; i++)
	{
		delete[]fg[i];
	}
	delete[]fg;
	printf("fgfcm complete\n\n");
	getchar();

	printf("start nlfcm\n");
	nlfcm(y, fg, width, height);

	statistic(fg, gt, width, height);
	for (i = 0; i<width; i++)
	{
		delete[]fg[i];
	}
	delete[]fg;
	printf("nlfcm complete\n\n");
	getchar();


	for (i = 0; i<width; i++)
	{
		delete[] x1[i];
		delete[] x2[i];
		delete[] y[i];
		delete[] gt[i];
	}
	delete[]x1;
	delete[]x2;
	delete[]y;
	delete[]gt;
}



void main()
{
	char s1[100], s2[100], s3[100];
	printf("请输入image1:");
	gets_s(s1);
	printf("请输入image2:");
	gets_s(s2);
	printf("请输入参考图:");
	gets_s(s3);

	//log(s1,s2,s3);
	//mean(s1,s2,s3);
	//wt(s1,s2,s3);
	int cl1, cl2;
	cl1 = clock();
	ct(s1, s2, s3);
	cl2 = clock();

	printf("%f\n", (double)(cl2 - cl1) / 1000.0);

	getchar();
	//comparafcm(s1,s2,s3);
}







	


