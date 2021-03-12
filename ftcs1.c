///////////////////
//////FTCS 1///////
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

	k = 0.4;
	tmax = 5000;

	for(i=0; i<L; i++)
	{
		f[i] = 0;
		g[i] = 0;
	}

	f[0] = 1;
	f[L-1] = 0;

	printf("set yrange [0:1]\n");
	printf("set key off\n");

	for(t=0; t < tmax+1; t++)
	{
		for(i=1; i<(L-1); i++)
		{
			g[i] = f[i] + k*(f[i-1] - 2.0*f[i] + f[i+1]);
		}

		for(i=1; i<(L-1); i++)
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
