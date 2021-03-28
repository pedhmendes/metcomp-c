/////////////////////////////
////DEFEITOS POR RADIAÇÃO////
////////////BTCS/////////////
////////PEDRO MENDES/////////
/////////////////////////////

//para rodar --> ./a.out | gnuplot

#include<stdio.h>
#include<math.h>
#include<stdlib.h>

#define L 350
#define s0 10
#define x0 100
#define sigma 10
#define S 0.0002

void tridiag_fixo(double a,double b,double c,double *d,double *u,int n,double dt);

int main()
{
	int tmax, i;
	double t, dt;
	double k;
	double *f_n, *f_n1;
	size_t size = L*sizeof(double);

	f_n = (double*)malloc(size);
	f_n1 = (double*)malloc(size);

	for(i=0; i<=L; i++)
	{
		f_n[i] = 0.;
		f_n1[i] = 0.;
	}	

	k = 0.4;	
	t = 0;
	dt = 0.25;
	tmax = 2000;

	printf("set key off\n");
	printf("set yrange[0:4000]\n");

	while(t<=tmax)
	{
		tridiag_fixo(-k, 1 + 2*k + S*dt, -k, f_n, f_n1, L,dt);
	
		for(i=0; i<L; i++)
		{
			f_n[i] = f_n1[i];// - s0*exp(-((pow((i-x0),2))/(2*pow(sigma,2))))*dt;
		}
		
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

void tridiag_fixo(double a,double b,double c,double *d,double *u,int n,double dt)
{
	double bet,*c_linha;
	int i, j;
	size_t size = n*sizeof(double);

	c_linha = (double*)malloc(size);

	for(i=0; i<=L; i++)
	{
		c_linha[i] = 0.;
	}	

	bet = b;

	for(j=1; j<n; j++)
	{
		d[j] = d[j] + s0*exp(-((pow((j-x0),2))/(2*pow(sigma,2))))*dt;
	}

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
