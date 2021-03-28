/////////////////////////////
////DEFEITOS POR RADIAÇÃO////
////////////FTCS/////////////
////////PEDRO MENDES/////////
/////////////////////////////

//para rodar --> ./a.out | gnuplot

#include<stdio.h>
#include<math.h>

#define L 350
#define s0 10
#define x0 100
#define sigma 10
#define S 0.0002

double s(double x)
{
	return s0*exp(-((pow((x-x0),2))/(2*pow(sigma,2))));
}

int main()
{
	int i, tmax;
	double f[L], g[L];
	double k, t, dt, Ddx;

	Ddx = 1.6;
	dt = 0.25;
	k = Ddx*dt;

	t=0;
	tmax = 2000;

	for(i=0; i<L; i++)
	{
		f[i] = 0;
		g[i] = 0;
	}

	printf("set yrange [0:4000]\n");
	printf("set key off\n");

	f[0] = 0;
	f[L-1] = 0;

	while(t <= tmax)
	{
		for(i=1; i<(L-1); i++)
		{
			g[i] = f[i] + k*(f[i-1] - 2.0*f[i] + f[i+1]) - (S*f[i] - s(i))*dt;
		}

		for(i=1; i<(L-1); i++)
		{
			f[i] = g[i];
		}

		printf("set title 'Tempo = %lf'\n", t);
		printf("plot \"-\" u 1:2 w lp\n");

		for(i=0; i<L; i++)
		{
			printf("%d	%lf\n", i, f[i]);
		}

		printf("e\n");

		t = t + dt;
	}

	printf("pause 5\n");

	return 0;
}
