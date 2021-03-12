///////////////////
//////FTCS 3///////
//////DIFUSAO//////
////PEDRO MENDES///
///////////////////

//para rodar --> ./a.out | gnuplot

#include<stdio.h>
#include<math.h>

#define L 100

int main()
{
	int i, t, tmax;
	double f[L], g[L];
	double k;
	int m1, m2;

	k = 0.4;
	tmax = 5000;

	for(i=0; i<L; i++)
	{
		f[i] = 0;
		g[i] = 0;
	}

	m1 = L/4;
	m2 = 3*(L/4);

	for(i=m1; i<=m2; i++)
	{
		f[i] = 1;
	}

	printf("set yrange [0:1]\n");
	printf("set key off\n");

	f[0] = 0;
	f[L-1] = 0;

	for(t=0; t < tmax+1; t++)
	{
		g[0] = f[0] + k*(f[L-1] - 2.0*f[0] + f[1]);
		g[L-1] = f[L-1] + k*(f[L-2] - 2.0*f[L-1] + f[0]);

		for(i=1; i<(L-1); i++)
		{
			g[i] = f[i] + k*(f[i-1] - 2.0*f[i] + f[i+1]);
		}

		for(i=0; i<L; i++)
		{
			f[i] = g[i];
		}

		printf("set title 'Tempo = %d'\n", t);
		printf("plot \"-\" u 1:2 w lp\n");

		for(i=0; i<L; i++)
		{
			printf("%d	%lf\n", i, f[i]);
		}

		printf("e\n");
	}

	printf("pause 5\n");

	return 0;
}
