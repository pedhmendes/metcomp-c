/********************************************************
 * 			TURING MODEL			*
 *******************************************************/ 		
/********************************************************
 * 	simulação do modelo de turing			*
 * 	modelo de reação-difusão			*
 * 							*
 * 	executar o programa inserindo   		*
 * 	os parametros na seguinte ordem 		*
 *	L a b c d Du Dv					*
 *							*
 *	exemplo de exec					*
 *	./a.out 100 1. -1. 2. -1.5 0.0001 0.0006	*
 *	./a.out L   a   b  c   d   Du	  Dv		*
 *							*
 *	se usar -DGNU tem que rodar 			*
 *	usando pipe pro gnuplot				*
 *	./a.out PARAMETROS | gnuplot			*
 *							*
 *	se usar -DDEBUG gera a mesma			*
 *	sequencia de numeros				*
 *******************************************************/ 

/********************************************************
 *			INCLUDES			*
 *******************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/********************************************************
 * 		      DEFINITIONS			*
 *******************************************************/ 
#define TMAX 5 

/********************************************************
 * 		   GLOBAL VARIABLES			*
 *******************************************************/ 
int L, L2;
double dh, dt;
double Du, Dv;
double a, b, c, d;
double h, k;

/********************************************************
 *	 	      FUNCTIONS				*
 *******************************************************/ 
void initialize(double *u, double *v, double *n_u, double *n_v, char *argv[ ]);
void update(double *u, double *v, double *n_u, double *n_v);
void gnuplot(double *vec);
double drand(double low, double high);

/********************************************************
 * 	 	    MAIN PROGRAM			*
 *******************************************************/ 
int main(int argc, char *argv[ ])
{
	double t;
	double *u, *v;
	double *n_u, *n_v;
	int count = 0;

	L = atoi(argv[1]);	
	L2 = L*L;
	
	size_t size = L*L*sizeof(double);
	
	u = (double*)malloc(size);
	v = (double*)malloc(size);
	n_u = (double*)malloc(size);
	n_v = (double*)malloc(size);

#ifdef DEBUG
	srand(time(0));	
#else
	srand(time(NULL));
#endif

	initialize(u,v,n_u,n_v,argv);	
	
	while(t < TMAX)
	{
		update(u,v,n_u,n_v);
#ifdef GNU
		if(count % 10 == 0)
		{
			printf("set title 'Tempo %.2lf'\n", t);
			gnuplot(v);
		}
#endif
		t = t + dt;
		count += 1;
	}

#ifdef GNU
	gnuplot(v);
	printf("pause 5\n");
#endif

	free(u);
	free(v);
	free(n_u);
	free(n_v);

	return 0;
}	

/********************************************************
 * 	   	      INITIALIZATION			*
 *******************************************************/ 
void initialize(double *u, double *v, double *n_u, double *n_v,char *argv[ ])
{
	for(int i=0; i<L2; i++)
	{
		u[i] = 1.0 + drand(-0.03, 0.03);	
		v[i] = 1.0 + drand(-0.03, 0.03);	
	
		n_u[i] = 0.0;
		n_v[i] = 0.0;
	}

	dh = 1.0/L;
	dt = 0.02;

	h = 1.;
	k = 1.;

	a = atof(argv[2]);
	b = atof(argv[3]);
	c = atof(argv[4]);
	d = atof(argv[5]);
	Du = atof(argv[6]);
	Dv = atof(argv[7]);

	return;
}

/********************************************************
 *	 		  UPDATE			*
 *******************************************************/ 
void update(double *u, double *v, double *n_u, double *n_v)
{
	int i;
	double uLap, vLap;

	for(i=0; i<L2; i++)
	{
		uLap = (u[(i-L+L2)%L2]+u[(i+1)%L + (i/L)*L]+u[(i+L)%L2]+u[(i-1+L)%L + (i/L)*L]-4*u[i])/(dh*dh);	
		vLap = (v[(i-L+L2)%L2]+v[(i+1)%L + (i/L)*L]+v[(i+L)%L2]+v[(i-1+L)%L + (i/L)*L]-4*v[i])/(dh*dh);

		n_u[i] = u[i] + (a*(u[i]-h) + b*(v[i]-k) +Du*uLap)*dt;	
		n_v[i] = v[i] + (c*(u[i]-h) + d*(v[i]-k) +Dv*vLap)*dt;	
	}
	
	for(i=0; i<L2; i++)
	{
		u[i] = n_u[i];
		v[i] = n_v[i];
	}	
	
	return;
}

/********************************************************
 * 			GNU VIEW			*
 *******************************************************/ 
void gnuplot(double *vec)
{        
	printf("set xrange [0:%d]\nset yrange [0:%d]\n", (L-1), (L-1));
	printf("unset xtics\nunset ytics\n");
	printf("unset colorbox\n");
	printf("set size square\n");
	printf("set key off\n");
	printf("plot \"-\" matrix w image\n");
	
	for(int i=0; i<L; i++)
	{
		for(int j=0; j<L; j++)
		{
			printf("%lf ", vec[i + j*L]);
		}
		printf("\n");
	}

	printf("\n\n\ne\n\n");

	return;
}

/********************************************************
 *		   RANDOM IN A RANGE            	*
 *******************************************************/
double drand(double low, double high)
{
	return ((double)rand()*(high-low))/(double)RAND_MAX + low;
}
