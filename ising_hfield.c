/*****************************************************************************
 *                             Ising HField                                  *
 *                            Pedro H Mendes                                 *
 ****************************************************************************/
/*****************************************************************************
 *                               INCLUDES                                    *
 ****************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/*****************************************************************************
 *                              DEFINITIONS                                  *
 ****************************************************************************/
#define 		L 			16
#define 		L2 			(L*L)
#define 		TRAN 			100000
#define 		TMAX			1000000
#define			HFIELD			0.005
#define			J			1.0

/*****************************************************************************
 *                           GLOBAL VARIABLES                                *
 ****************************************************************************/
int site[L][L], neigh[L][L];
double M, ET;

/*****************************************************************************
 *                              FUNCTIONS                                    *
 ****************************************************************************/
void initialize(void);
void routine(double TEMP);
void gnu_view(void);

/*****************************************************************************
 *                             MAIN PROGRAM                                  *
 ****************************************************************************/
int main(int argc, char *argv[])
{
	int mcs;
	double TEMP;
	char Arq1[100];
	FILE *arq1;

	TEMP = atof(argv[1]);

	srand(time(NULL));
	initialize();

	for(mcs=0; mcs<TRAN; mcs++)
	{
		routine(TEMP);
	}

	sprintf(Arq1, "series_L%d_T%.3lf.dsf", L, TEMP);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#MCS\tM\tET\n");
	
	for(mcs=0; mcs<TMAX; mcs++)
	{
		routine(TEMP);
		fprintf(arq1, "%d\t%lf\t%lf\n", mcs, M, ET);
#ifdef GNU
		if(mcs%10 == 0)
		{
			printf("set title 'MCS = %d'\n", mcs);
			gnu_view();
		}
#endif
	}

	fclose(arq1);

	return 0;
}

/*****************************************************************************
 *                            INITIALIZATION                                 *
 ****************************************************************************/
void initialize(void)
{
	int i, j;
	int E1, E2;

	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			site[i][j] = 2*(rand()%2) - 1;
			neigh[i][j] = 0;
		}
	}

	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			neigh[i][j] += site[(i == 0 ? (L-1) : i - 1)][j];
			neigh[i][j] += site[(i == (L-1) ? 0 : i + 1)][j];
			neigh[i][j] += site[i][(j == 0 ? (L-1) : j - 1)];
			neigh[i][j] += site[i][(j == (L-1) ? 0 : j + 1)];
		}
	}

	ET = 0;
	E1 = 0;
	E2 = 0;
	M = 0;

	for(i=0; i<L; i++)
	{
		for(j=0; j<L; j++)
		{
			E1 = E1 + site[i][j]*neigh[i][j];
			E2 = E2 + site[i][j];
			M = M + site[i][j];
		}
	}

	E1 = E1*0.5;

	ET = -E1-(E2*HFIELD);		

	return;
}

/*****************************************************************************
 *                            MONTE CARLO ROUTINE                            *
 ****************************************************************************/
void routine(double TEMP)
{
	int i, j, t;
	double dE;
	double r, prob;

	for(t=0; t<L2; t++)
	{
		i = rand()%L;
		j = rand()%L;
	
		dE = 0;
		dE = 2.0*site[i][j]*(J*neigh[i][j] + HFIELD);    

		if(dE <= 0)
		{
			site[i][j] *= -1;
			ET = ET + dE;
			M = M + 2*site[i][j];

			neigh[(i == 0 ? (L-1) : i - 1)][j] += 2*site[i][j];
			neigh[(i == (L-1) ? 0 : i + 1)][j] += 2*site[i][j];
			neigh[i][(j == 0 ? (L-1) : j - 1)] += 2*site[i][j];
			neigh[i][(j == (L-1) ? 0 : j + 1)] += 2*site[i][j];
		}
		else
		{
			prob = exp(-dE/TEMP);
			r = 1.0*rand()/RAND_MAX;
		
			if(r < prob)
			{
				site[i][j] *= -1;
				ET = ET + dE;
				M = M + 2*site[i][j];

				neigh[(i == 0 ? (L-1) : i - 1)][j] += 2*site[i][j];
				neigh[(i == (L-1) ? 0 : i + 1)][j] += 2*site[i][j];
				neigh[i][(j == 0 ? (L-1) : j - 1)] += 2*site[i][j];
				neigh[i][(j == (L-1) ? 0 : j + 1)] += 2*site[i][j];
			}
		}
	}

	return;
}

/*****************************************************************************
 *                             GNUPLOT VIEW                                  *
 ****************************************************************************/
void gnu_view(void)
{
	int a, b;

	printf("set xrange [0:%d]\nset yrange [0:%d]\n", (L-1), (L-1));
	printf("set cbrange [0:1]\n");
	printf("unset xtics\nunset ytics\n");
	printf("set palette defined (-1 'black', 1 'white')\n");   
	printf("unset colorbox\n");
	printf("set size square\n");
	printf("set key off\n");
	printf("plot \"-\" matrix w image\n");

	for(a=0; a<L; a++)
	{
		for(b=0; b<L; b++)
		{
			printf("%d ", site[a][b]);
		}
		printf("\n");
	}

	printf("\ne\n\n");
	printf("pause 0.1\n\n");

	return;
}
