/*#include <ostream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
*/
#include <stdlib.h>
#include <stdio.h>
#include <fstream.h>
#include <malloc.h>
#include <time.h>
#include <math.h>


#define debug 1
#define maxn  10000    //for poisson  self-similar zipf

int cnum=1;

float Uniform();
float Gauss(double mean, double deviat);  //mean=0.5, deviat=1/6
float ExpNeg(float lamder);  //lamder=0.2
float Poisson();  //lamder=3330
float SelfSim(double h); // 0<h<1
float Zipf(long n, double theta);

float Cluster(int clusternum, float *mean, float *deviate)
{
	float *temp=new float[clusternum];
//	float *mean=new float[clusternum];
//	float *deviate=new float[clusternum];

	for (int i=0; i<clusternum; i++)
	{
//		srand((unsigned)clock());
//		mean[i]=Uniform();
//		srand((unsigned)clock());
//		deviate[i]=Uniform();
		temp[i]=Gauss(mean[i], deviate[i]);
	}
//	srand((unsigned)clock());
	i=rand()%clusternum;

	return temp[i];
}
	
/*Uniform() */
float Uniform()
{
	return ((double)rand()/(double)RAND_MAX);
}
/* Gauss() */
/*
Gauss(DistrAT * distr)
{
    int				i;
    double			ans = 0;
    for (i = 0; i < 12; i++) {
        ans += distr->iRand /(double) distr->Prim1 - 0.5;
        distr->iRand = MulM(distr->iRand, distr->Root, distr->Prim);
    }
    return (__int64) (distr->Mean + distr->Deviat * ans);
}
*/
/* Gauss() */
float Gauss(double mean, double deviat)  //mean=0.5, deviat=1/6
{
    int				i;
    double			ans = 0;
    for (i = 0; i < 12; i++) {
        ans += (double)rand() /(double)RAND_MAX  - 0.5;
    }
    return (mean + deviat * ans);
}

/* ExpNeg()*/
/*
__int64			ExpNeg(DistrAT * distr)
{
    __int64			ret;
    ret = (__int64) (- distr->Deviat * log(distr->iRand / (double) distr->Prim1));
    distr->iRand = MulM(distr->iRand, distr->Root, distr->Prim);
    return ret;
}*/

/* ExpNeg()*/
float ExpNeg(double lamder)  //lamder=1/6
{
    return (-(log((double)rand())-log((double)RAND_MAX))*lamder);
}

/* Poisson() */
/*
__int64			Poisson(DistrAT * distr)
{
    __int64			i = 0;
    double			c, p = 1;
    c = pow(2.718282, - distr->Deviat);
    if (c == 0)
        return 0;
    while (p >= c) {
        p = p * (distr->iRand / (double) distr->Prim1);
        i++;
		distr->iRand = MulM(distr->iRand, distr->Root, distr->Prim);
    }
    return (i - 1);
}*/
float Poisson()  
{
    long	i = 0;
    double	c, p = 1;
	
	double lamder=2000.0/3.0;

    c = pow(2.718282, -(double)lamder);
    if (c == 0)
        return 0;
    while (p >= c) {
        p = p * ((double)rand() / RAND_MAX);
        i++;
    }
    return ((double)(i - 1)/(double)2000);
}

/* SelfSim() */
/*
__int64			SelfSim(DistrAT * distr)
{
    __int64			ret;
    ret = (__int64) (1 + distr->iRand * pow(distr->iRand / (double) distr->Prim1,
        log(distr->Deviat) / log(1 - distr->Deviat)));
    distr->iRand = MulM(distr->iRand, distr->Root, distr->Prim);
    return ret;
}*/
float SelfSim(double h) // 
{
    double ret;
    ret = (1 + (int)(maxn * pow((double)rand()/(double)RAND_MAX, log(h) / log(1 - h))))/(double)maxn;
    
    return ret;
}
/* zeta() */
/*
double			zeta(__int64 n, double theta)
{
    int				i;
    double			ret = 0.0;
    for (i = 1; i < n; i++)
        ret += pow(i, -theta);
    return ret;
}*/
double zeta(int n, double theta)
{
    int				i;
    double			ret = 0.0;
    for (i = 1; i < n; i++)
        ret += pow(i, -theta);
    return ret;
}

/* Zipfr() */
/*
	__int64			Zipf(DistrAT * distr)
{
    __int64			r, ret;
    double			alpha, zetan, eta, uz, u;
    r = distr->iRand;
    distr->iRand = MulM(distr->iRand, distr->Root, distr->Prim);
    alpha = 1 / (1 - distr->Deviat);
    u = r / (double) distr->Prim1;
    zetan = zeta(r, distr->Deviat);
    uz = u * zetan;
    if (uz < 1)
        return 1;
    if (uz < 1 + pow(0.5, distr->Deviat))
        return 2;
    eta = (1 - pow(2.0 / r, 1 - distr->Deviat)) * (1 - zeta(2, distr->Deviat) / zetan);
    ret = (__int64) (1 + r * pow(eta * u - eta + 1, alpha));
    return ret;
}*/

float Zipf(double theta)
{
	double u=(double)rand()/(double)RAND_MAX;
    double zetan = zeta(maxn, theta);
    double uz = u * zetan;

    if (uz < 1)
        return (1.0/(double)maxn);
    
	if (uz < 1 + pow(0.5, theta))
        return (2.0/(double)maxn);

	double alpha = 1 / (1 - theta);
    double eta = (1 - pow(2.0 /(double) maxn, 1- theta))*(1 - zeta(2, theta) / zetan);
    double ret =(double) ((int)(1 + maxn * pow(eta * u - eta + 1, alpha)))/(double)maxn;
    return ret;
}

float gendata(int distribution, int i, float *mean, float *deviate)
{
	float iVal;

        do {
            switch (distribution) {
            case 0: //ordered
				iVal=i; 
                break;
            case 1: //DISTR_Uniform:
				iVal=Uniform(); 
				break;
            case 2: //DISTR_GAUS: Gauss(float mean, float deviat)  //mean=0.5, deviat=1/6
//				double deviat=1.0/6.0;
                iVal = Gauss(0.5, 1.0/6.0);
                break;
            case 3: //DISTR_EXPN:ExpNeg(float lamder);  //lamder=1/6
                iVal = ExpNeg(1.0/6.0);
                break;
            case 4: //DISTR_POIS: Poisson(float lamder);  //lamder=3333
                iVal = Poisson();
                break;
            case 5: //DISTR_SSIM:SelfSim(long n, double h); // 0<h<1
                iVal = SelfSim(0.1);
                break;
            case 6: //DISTR_ZIPF:Zipf(long n, double theta);
                iVal = Zipf(0.5);
                break;
			case 7:
				iVal=Cluster(cnum, mean, deviate);
				break;
            }
		}
			while (iVal > 1 || iVal < 0);
//#ifdef debug
//		printf("value: %f\n",iVal);
//#endif
        return(iVal);
}


int	main(int argc, char *argv[])
{
    int	i;
	float *mean;
	float *deviate;

    if (argc < 5) {
        printf("usage: gdt <dimension> <number> <distribution > <output file>\n");
		printf("Distribution: <0 ordered> <1 uniform> <2 Normal> <3 Exp> \n");
		printf("              <4 Poisson> <5 selfsim> <6 Zipf> <7 cluster>\n");
		printf("Usage for generate cluster: gdt <dim> <num> <dis> <file> <cnum>\n");
        return 0;
    }

	int dimension=2;//atoi(argv[1]);
	int number=5000;//atoi(argv[2]);
	int distribution=1;//atoi(argv[3]);


	if (distribution==7){
		if (argc!=6){
		printf("Usage for generate cluster: gdt <dim> <num> <dis> <file> <cnum>\n");
        return 0;}
		else 
		{
			cnum=atoi(argv[5]);
			mean=new float[cnum];
			deviate=new float[cnum];
			//srand((unsigned int)clock());
			for (int i=0; i<cnum;i++)
				{
					mean[i]=Uniform();
					deviate[i]=0.05;
					//deviate[i]=Uniform();
				}
			for (i=0; i<cnum;i++)
			{ 
				cout << "mean["<<i<<"]"<< mean[i] << "\t";
			}
			cout << endl;
			getchar();
		}
	}
	else { 
		if (argc!=5) {
        printf("usage: gdt <dimension> <number> <distribution > <output file>\n");
		printf("Distribution: <0 ordered> <1 uniform> <2 Normal> <3 Exp> \n");
		printf("              <4 Poisson> <5 selfsim> <6 Zipf> <7 cluster>\n");
		printf("Usage for generate cluster: gendata <dim> <num> <dis> <file> <cnum>\n");
        return 0;
		}
		else  {
			float *mean=NULL;
			float *deviate=NULL;
		}
	}


	float *buffer=new float [number*dimension];

	if (buffer==NULL)
	{
		printf("unable to allocate memory!\n");
		return (0);
	}

	
/*	
	if ((buffer=malloc(number*dimension*sizeof(float)))==NULL) 
	{
		printf("unable to allocate memory!\n");
		return (0);
	}
*/
    for (i = 0; i < dimension; i++) {
		//srand( (unsigned)rand() );
		for (int j=0, x=i; j<number; j++, x+=dimension) 
		{
			//cout << "i:" <<i << "  j:"<<j<<endl;
			buffer[x]=gendata(distribution,j, mean, deviate);
		}
    }

	cout << "writing to " << argv[4]<<endl;
	FILE *fp=fopen("C:/LBSQ/Datasets/uni5000.txt","wt");

	if (fp==NULL) {
		cout << "cannot creat file "<<argv[4]<<endl;
		return(0);
	}
	
	float *floatp=(float*)buffer;
	int j;
	for (i=0; i<number; i++)
	{
		for (j=0;j<dimension;j++)
		{
			fprintf(fp, "%f\t", *floatp);
			floatp++;
		}
		fprintf(fp, "\n");
	}
	
	cout << "write "<<i<<" number"<<endl;
	fclose(fp);
	delete buffer;

	cout <<"end!"<<endl;

	return (0);

}

