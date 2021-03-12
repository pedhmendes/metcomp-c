///////////////////////
////CRANK-NICOLSON/////
///////DIFUSAO 2///////
/////PEDRO MENDES//////
///////////////////////

//para rodar --> ./a.out | gnuplot

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "pointers.h"

#define L 100

void tridiag_fixo(double a,double b,double c,double *d,double *u,int n);

int main()
{
	int tmax, i;
	double m1, m2;
	double t, dt, k;
	double *f_n, *g_n;
	size_t size = (L*sizeof(double)) + 1.;

	f_n = (double*)malloc(size);
	g_n = (double*)malloc(size);

	for(i=0; i<=L; i++)
	{
		f_n[i] = 0.;
		g_n[i] = 0.;
	}	

	k = 0.4;
	t = 0;
	dt = 0.4;

	m1 = L/4;
	m2 = 3*(L/4);
	tmax = 2000;

	for(i=m1; i<=m2; i++)
	{
		f_n[i] = 1;
	}

	printf("set key off\n");
	printf("set yrange[0:1]\n");

	while(t<=tmax)
	{
		for(i=1; i<(L-1); i++)
		{
			g_n[i] = k*f_n[i-1]+(2-2*k)*f_n[i]+k*f_n[i+1];
		}

		tridiag_fixo(-k, 2 + 2*k, -k, g_n, f_n, L+1);
	
		printf("set title 't = %lf'\n", t);
		printf("plot \"-\" u 1:2 w lp\n");

		for(i=0; i<=L; i++)
		{
			printf("%d	%lf\n", i, f_n[i]);
		}
		
		printf("e\n");

		t = t + dt;
	}

	printf("pause 5\n");

	return 0;
}

void tridiag_fixo(double a,double b,double c,double *d,double *u,int n)
{
	double bet,*c_linha;
	int j;
	
	c_linha = create_double_pointer(n);
	bet = b;

	d[1] = d[1] - a*d[0];
	u[1] = d[1]/bet;
	d[n-2] = d[n-2] - c*d[n-1];
	c_linha[1] = c/b;
	
	for(j=2;j<n-1;j++)
	{
    		bet = b - a*c_linha[j-1];
    		c_linha[j] = c/bet;
    		u[j] = (d[j]-a*u[j-1])/bet;
	}

	for(j=(n-2);j>0;j--)
	{
    		u[j] = u[j] - c_linha[j]*u[j+1];
  	}

	free(c_linha);
}
