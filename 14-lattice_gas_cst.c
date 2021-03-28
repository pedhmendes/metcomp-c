/****************************************************************************
 *    			     Lattice Gas Cst				    *
 * 			      Pedro H Mendes				    *
 ***************************************************************************/

/****************************************************************************
 * 			    	 INCLUDES				    *
 ***************************************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/****************************************************************************
 * 			        DEFINITIONS				    *
 ***************************************************************************/
#define 			L				100
#define				L2				(L*L)
#define 			RHO				0.25
#define 			EPSILON				(1.0)
#define				TRAN				1000
#define 			TMAX				10000

/****************************************************************************
 * 			     GLOBAL VARIABLES				    *
 ***************************************************************************/
int site[L2], neigh[L2][5];
int ET, dE;
int density;
double boltz[4];

/****************************************************************************
 * 			        FUNCTIONS				    *
 ***************************************************************************/
void initialize(double TEMP);
void monte_carlo_routine(void);
void update_neighbour(int k);
void density_function(void);
void gnu_plot(void);

/****************************************************************************
 * 			       MAIN PROGRAM				    *
 ***************************************************************************/
int main(int argc, char *argv[])
{
	int mcs;
	double TEMP;
	char Arq1[100];
	FILE *arq1;

	TEMP = atof(argv[1]);

	srand(time(NULL));

	initialize(TEMP);

	for(mcs=0; mcs<TRAN; mcs++)
	{
		monte_carlo_routine();
		if(mcs%10 == 0) 
	//	{
	//		printf("set title 'Tempo = %d MCS\n", mcs);
			gnu_plot();
	//	}	
	}

	sprintf(Arq1, "series_T%.3lfL%d.dsf", TEMP, L);
	arq1 = fopen(Arq1, "w");
	fprintf(arq1, "#MCS\tET\tdensity\n");

	for(mcs=0; mcs<TMAX; mcs++)
	{
		monte_carlo_routine();
		density_function();
		fprintf(arq1, "%d\t%d\t%d\n", mcs, ET, density);
#ifdef GNU
		if(mcs%10 == 0) 
	//	{
	//		printf("set title 'Tempo = %d MCS\n", mcs);
			gnu_plot();
	//	}
#endif

	}
	
	fclose(arq1);

	return 0;
}

/****************************************************************************
 * 			      INITIALIZATION				    *
 ***************************************************************************/
void initialize(double TEMP)
{
	int i, j;

	boltz[1] = exp(-1.0/TEMP);
	boltz[2] = exp(-2.0/TEMP);
	boltz[3] = exp(-3.0/TEMP);

	for(i=0; i<L2; i++)
	{
		site[i] = 0;
		
		for(j=0; j<5; j++)
		{	
			neigh[i][j] = 0;
		}
	}

	for(i=0; i<(RHO*L2); i++)  
	{
		j = rand()%L2;	
		site[j] = 1;
	}
	
	ET = 0;
	density = 0;

	for(i=0; i<L2; i++)
	{
		neigh[i][0] = (i-L+L2)%L2;
		neigh[i][1] = (i+1)%L + (i/L)*L;
		neigh[i][2] = (i+L)%L2;
		neigh[i][3] = (i-1+L)%L + (i/L)*L;

		for(j=0; j<4; j++)
		{
			neigh[i][4] += site[neigh[i][j]];
		}

		ET += site[i]*neigh[i][4];

		density += site[i];
	}

	ET = -ET/2;

	return;
}


/****************************************************************************
 * 			      MONTE CARLO ROUTINE			    *
 ***************************************************************************/
void monte_carlo_routine(void)
{
	int i, j, k, t;
	int s_i, s_j;
	int sum_i, sum_j;
	int site_aux;
	int Ei, Ef;
	double r, flip_prob;

	for(t=0; t<L2; t++)
	{
		i = rand()%L2;
		k = rand()%4;
		j = neigh[i][k];

		s_i = site[i];
		s_j = site[j];
		
		if(s_i != s_j)
		{
			sum_i = neigh[i][4];
			sum_j = neigh[j][4];
	
			Ei = -sum_i*s_i - sum_j*s_j;

			sum_i = sum_i + s_i - s_j;
			sum_j = sum_j + s_j - s_i; 
	
			Ef = -sum_i*s_j - sum_j*s_i;
	
			dE = 0;
			dE = Ef-Ei;

			if(dE <= 0)
			{
				site_aux = site[i];
				site[i] = site[j];
				site[j] = site_aux;
	
				ET += dE;
	
				update_neighbour(i);
				update_neighbour(j);
			}
			else
			{
				flip_prob = boltz[dE];
				r = rand()/(double)RAND_MAX;
	
				if(r < flip_prob)
				{
					site_aux = site[i];
					site[i] = site[j];
					site[j] = site_aux;
						
					ET += dE;
	
					update_neighbour(i);
					update_neighbour(j);
				}
			}

		}
	}

	return;
}


/****************************************************************************
 * 			     UPDATE NEIGHBOURS			  	    *
 ***************************************************************************/
void update_neighbour(int k)
{
	int i;
	
	for(i=0; i<4; i++)
	{
		neigh[neigh[k][i]][4] += 2*site[k] -1;
	}

	return;
}

/****************************************************************************
 * 			     DENSITY FUNCTION			  	    *
 ***************************************************************************/
void density_function(void)
{
	int i;

	density = 0;

	for(i=0; i<L2; i++)
	{
		density += site[i];
	}

	return;
}

/****************************************************************************
 * 		    	       GNUPLOT VIEW		 	 	    *
 ***************************************************************************/
void gnu_plot(void)
{
	int ii,jj;

	printf("set xrange [0:%d]\nset yrange [0:%d]\n", (L-1), (L-1));
	printf("set cbrange [0:4]\n");
	printf("unset xtics\nunset ytics\n");
	printf("set palette defined (0 'black', 1 'dark-violet')\n");
	printf("unset colorbox\n");
	printf("set size square\n");
	printf("set key off\n");
	printf("plot \"-\" matrix w image\n");

	for(jj=0;jj<L;jj++)
	{
		for(ii=0;ii<L;ii++)
		{
			printf("%d ", site[ii + jj*L]);
		}
		printf("\n");
	}

//	printf("\ne\n pause 0.1\n");
//	printf("\ne\n pause 0.000001\n");
	printf("\ne\n\n");


	return;
}
