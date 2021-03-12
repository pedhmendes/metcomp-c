/////////////////////////////////////////////////////////////////////////
////////////////////////DILEMA DO PRISIONEIRO////////////////////////////
///////////////////////////PEDRO MENDES//////////////////////////////////
/////////////////////////////////////////////////////////////////////////

/************************************************************************
 *	 	Requer a bilioteca cientifica GSL			*
 *									*
 * Para compilar a versao mais simples do programa use:                 *
 *->gcc -Wall prisoner_dilemma.c -lgsl -lgslcblas -lm -static           *
 *  									*
 * A execução requer que insira o valor de custo beneficio na linha de 	*
 * de comando da forma:							*
 *->./a.out r								*
*									*
 * Para ver um time lapse do sistema ao longo da simulação compile:	*
 *->gcc -Wall -DLAPSE prisoner_dilemma.c -lgsl -lgslcblas -lm -static	*
 *									*
 * Se compilar com a opcao de timelapse do sistema deve rodar com um 	*
 * pipe para gnuplot da forma:						*
 *-> ./a.out r | gnuplot						*
 *									*
 * Dependendo dos tamanhos de L, transient e tmax pode levar um tempo	*
 * para executar, cuidado						*
 ***********************************************************************/

/*******INCLUDES*******/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_rng.h>

/*******DEFINICOES*******/
#define L 1000
#define L2 (L*L)
#define transient 100000 		// 1e5
#define tmax (1000000-transient)	// 1e6
#define NV 4

/*******VARIAVEIS GLOBAIS*******/
int x[L][L];
int n[L][L];
int hist_coop[L2];
int coop;

/*******FUNCOES*******/
double f_prob(double payx, double payy);
void coop_func(void);
void initialize(gsl_rng *rr);
void gnu_plot(double time);
void snapshot(FILE *file);
void hist_plot(FILE *file);

/*******MAIN*******/
int main(int argc, char *argv[])
{
	double T, R, P, S;
	double r;
	double payC, payD, payX, payY;
	double rnd, wprob;
	int mcs, tt, i, ii=0, j, jj=0;
	int c;
	double coop_media;
	char Serie_T[1000], Snap_File[1000], H_Coop[1000];
	FILE *serie_t, *snap_file, *h_coop;
	unsigned long int seed;

/*******SEMENTE DO GERADOR*******/
	srand(time(NULL));
	seed = rand();

/*******GERADOR DE NUMEROS ALEATORIOS*******/
	gsl_rng *rr; 
	rr = gsl_rng_alloc(gsl_rng_mt19937); 
	gsl_rng_set(rr,seed); 


/*******RELACOES DE CUSTO-BENEFICIO*******/
	r = atof(argv[1]);
	T = (1+r);
	R = 1.0;
	P = 0.0;
	S = (-1*r);

/*******INICIALIZACAO DO SISTEMA*******/
	initialize(rr); 

/*******COMECA SIMULACAO COM TEMPO TRANSIENTE*******/	
	for(mcs=0; mcs<transient; mcs++)
	{

/*******PASSO DE MONTE CARLO POR SITIO*******/
		for(tt=0; tt<L2; tt++) 
		{

/*******SELECIONA SITIO E CALCULA PAYOFF*******/
			i = gsl_rng_uniform_int(rr,L); 
			j = gsl_rng_uniform_int(rr,L); 
	
			payC = n[i][j]*R + (NV-n[i][j])*S; 
			payD = n[i][j]*T + (NV-n[i][j])*P; 
			payX = x[i][j]*payC + (1-x[i][j])*payD; 

/*******SELECIONA VIZINHO E CALCULA O PAYOFF*******/
			c = gsl_rng_uniform_int(rr,4);

			switch (c)
			{ 
				case 0:
					ii = (i==0 ? (L-1) : i-1);jj = j;
				break;
				case 1:
					ii = (i==(L-1) ? 0 : i+1);jj = j;
				break;
				case 2:	
					ii = i;jj = (j==0 ? (L-1) : j-1);
				break;
				case 3:
					ii = i;jj = (j==(L-1) ? 0 : j+1);
				break;
			}

			payC = n[ii][jj]*R + (NV-n[ii][jj])*S; 
			payD = n[ii][jj]*T + (NV-n[ii][jj])*P; 
			payY = x[ii][jj]*payC + (1-x[ii][jj])*payD; 
			
/*******CALCULA PROBABILIDADE DE MUDAR DE ESTRATEGIA*******/
			wprob = f_prob(payX, payY); 
			rnd = gsl_rng_uniform(rr); 
			
			if(wprob > rnd)
			{
/*******SE MUDAR COPIA A ESTRATEGIA DO VIZINHO*******/
				if(x[i][j] != x[ii][jj])
				{	
					x[i][j] = x[ii][jj]; 
					n[(i==0 ? (L-1) : i-1)][j] += 2*(x[i][j]) - 1;
					n[(i==(L-1) ? 0 : i+1)][j] += 2*(x[i][j]) - 1;
					n[i][(j==0 ? (L-1) : j-1)] += 2*(x[i][j]) - 1;
					n[i][(j==(L-1) ? 0 : j+1)] += 2*(x[i][j]) - 1;
				}
			}
		}

/*******TIMELAPSE*******/
#ifdef LAPSE
		printf("set title '%d MCS'\n", mcs);
		gnu_plot(0.000005);
#endif
	}

/*******AGORA CALCULA E SALVA DADOS*******/
	sprintf(Serie_T, "serie_r%lf.dat", r);
	serie_t = fopen(Serie_T,"w"); 

	for(mcs=0; mcs<tmax; mcs++)
	{
/*******PASSO DE MONTE CARLO POR SITIO*******/
		for(tt=0; tt<L2; tt++) 
		{
/*******SELECIONA SITIO E CALCULA PAYOFF*******/
			i = gsl_rng_uniform_int(rr,L); 
			j = gsl_rng_uniform_int(rr,L); 
	
			payC = n[i][j]*R + (NV-n[i][j])*S; 
			payD = n[i][j]*T + (NV-n[i][j])*P; 
			payX = x[i][j]*payC + (1-x[i][j])*payD; 

/*******SELECIONA VIZINHO E CALCULA O PAYOFF*******/
			c = gsl_rng_uniform_int(rr,4);

			switch (c)
			{ 
				case 0:
					ii = (i==0 ? (L-1) : i-1);jj = j;
				break;
				case 1:
					ii = (i==(L-1) ? 0 : i+1);jj = j;
				break;
				case 2:	
					ii = i;jj = (j==0 ? (L-1) : j-1);
				break;
				case 3:
					ii = i;jj = (j==(L-1) ? 0 : j+1);
				break;
			}

			payC = n[ii][jj]*R + (NV-n[ii][jj])*S; 
			payD = n[ii][jj]*T + (NV-n[ii][jj])*P; 
			payY = x[ii][jj]*payC + (1-x[ii][jj])*payD; 
			
/*******CALCULA PROBABILIDADE DE MUDAR DE ESTRATEGIA*******/
			wprob = f_prob(payX, payY); 
			rnd = gsl_rng_uniform(rr); 
			
			if(wprob > rnd)
			{
/*******SE MUDAR COPIA A ESTRATEGIA DO VIZINHO*******/
				if(x[i][j] != x[ii][jj])
				{	
					x[i][j] = x[ii][jj]; 
					n[(i==0 ? (L-1) : i-1)][j] += 2*(x[i][j]) - 1;
					n[(i==(L-1) ? 0 : i+1)][j] += 2*(x[i][j]) - 1;
					n[i][(j==0 ? (L-1) : j-1)] += 2*(x[i][j]) - 1;
					n[i][(j==(L-1) ? 0 : j+1)] += 2*(x[i][j]) - 1;
				}
			}
		}

/*******CALCULA A QUANTIDADE DE COOPERADORES E ATUALIZA HISTOGRAMA*******/
		coop_func();  
		hist_coop[coop] += 1;
		fprintf(serie_t,"%d %d %d\n", mcs, coop, (L2-coop));
		
/*******TIMELAPSE*******/
#ifdef LAPSE
		printf("set title '%d MCS'\n", mcs);
		gnu_plot(0.000005);
#endif
	}
	fclose(serie_t);
	
/*******FAZ A MEDIA DO HISTOGRAMA PRA ACHAR A MEDIA DE COOP*******/
	coop_media = 0.;

	for(i=0; i<L2; i++)
	{
		coop_media += i*hist_coop[i];
	}

	coop_media = (1.0*coop_media/tmax);

	printf("%lf\t%lf\t%lf\n", r, coop_media, (L2-coop_media));

/*******IMPRIME HISTOGRAMA DE COOPERADORES*******/
	sprintf(H_Coop, "hist_c_r%lf.dat", r);
	h_coop = fopen(H_Coop,"w"); 
	hist_plot(h_coop);
	fclose(h_coop);

/*******SNAPSHOT FINAL DO SISTEMA*******/
	sprintf(Snap_File, "fcond_r%lf.gnu", r);
	snap_file = fopen(Snap_File, "w"); 
	snapshot(snap_file);	
	fclose(snap_file);

/*******LIMPA MEMORIA DO GERADOR*******/
	gsl_rng_free(rr); 

	return 0;
}

/*******EQUACAO DE PROBABILIDADE*******/
double f_prob(double payx, double payy)
{
	double k=0.1;

	return 1/(1+exp((payx-payy)/k));
}

/*******CALCULA QUANTIDADE DE COOPERADORES POR MCS*******/
void coop_func(void)
{
	int i, j;

	coop=0;

	for(i=0; i<L; i++) 
	{
		for(j=0; j<L; j++)
		{
			if(x[i][j] == 1)
			{			
				coop = coop + 1;
			}
		}
	}

	return;
}

/*******INICIALIZA O SISTEMA*******/
void initialize(gsl_rng *rr)
{
	int i ,j;

/*******CONDICAO ALEATORIA DE ESTRATEGIAS*******/
	for(i=0; i<L; i++) 
	{
		for(j=0; j<L; j++)
		{
			x[i][j] = gsl_rng_uniform_int(rr,2);
		}
	}

/*******CRIA MATRIZ DE VIZINHOS*******/
	for(i=0; i<L; i++) 
	{
		for(j=0; j<L; j++)
		{
			n[i][j] = 0;
			n[i][j] += x[(i == 0 ? (L-1) : i - 1)][j];
			n[i][j] += x[(i == (L-1) ? 0 : i + 1)][j];
			n[i][j] += x[i][(j == 0 ? (L-1) : j - 1)];
			n[i][j] += x[i][(j == (L-1) ? 0 : j + 1)];
		}
	}

/*******LIMPA ARRAY DO HISTOGRAMA*******/
	for(i=0; i<L2; i++)
	{
		hist_coop[i] = 0;
	}

	return;
}

/*******FUNCAO PRA TIMELAPSE DO SISTEMA*******/
void gnu_plot(double time)
{
	int ii,jj;

	printf("set xrange [0:%d]\nset yrange [0:%d]\n", (L-1), (L-1));
	printf("set cbrange [0:1]\n");
	printf("unset xtics\nunset ytics\n");
	printf("set palette grey negative\n");
	printf("unset colorbox\n");
	printf("set size square\n");
	printf("set key off\n");
	printf("plot \"-\" matrix w image\n");
	
	for(jj=0;jj<L;jj++)
	{
		for(ii=0;ii<L;ii++)
		{
			printf("%d ", x[ii][jj]);
		}
		printf("\n");
	}
	printf("\ne\npause %lf\n\n\n\n", time);
	
	return;
}

/*******SNAPSHOT FINAL DO SISTEMA*******/
void snapshot(FILE *file)
{
	int ii, jj;
	
	fprintf(file, "set xrange [0:%d]\nset yrange [0:%d]\n", (L-1), (L-1));
	fprintf(file, "set cbrange [0:1]\n");
	fprintf(file, "unset xtics\nunset ytics\n");
	fprintf(file, "set palette grey negative\n");
	fprintf(file, "unset colorbox\n");
	fprintf(file, "set size square\n");
	fprintf(file, "set key off\n");
	fprintf(file, "plot \"-\" matrix w image\n");
	
	for(jj=0;jj<L;jj++)
	{	
		for(ii=0;ii<L;ii++)
		{
			fprintf(file, "%d ",x[ii][jj]);
		}
		fprintf(file, "\n");
	}

	fprintf(file, "e\n");

	return;
}

/*******FUNCAO PRA PLOTAR HISTOGRAMA*******/
void hist_plot(FILE *file)
{	
	int i;
	int i_min, i_max;
	
	fprintf(file, "#i\thist_coop[i]\n");
	
	i_min = 0;
	i_max = L2-1;

/*******MENOR VALOR DE RELEVANCIA*******/
	for(i=0; i<L2; i++)
	{
		if(hist_coop[i] != 0)
		{
			i_min = i;
			break;
		}
	}

/*******MAIOR VALOR DE RELEVANCIA*******/
	for(i=(L2-1); i>=0; i--)
	{
		if(hist_coop[i] != 0)
		{
			i_max = i;
			break;
		}
	}

/*******IMPRIME DENTRO DO INTERVALO*******/
	for(i=i_min; i<=i_max; i++)
	{
		fprintf(file, "%d\t%d\n", i, hist_coop[i]);
	}
	
	return;
}
