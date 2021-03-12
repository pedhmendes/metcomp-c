///////////////////////////
////////PEDRO MENDES///////
///////////////////////////

// compile	-> gcc -Wall ho_euler_gsl.c -lgsl -lgslcblas -lm
// run		-> ./a.out > serie.dsf
	
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define 	TMAX		100 
#define 	dt		0.001
#define 	ALPHA		1
#define		BETA		1	
#define		GAMMA		1

int main()
{
	double t, x, v;
	double n_gauss, sigma;
	unsigned long int seed;

	srand(time(NULL));
	seed = rand();

	gsl_rng *mt;
	mt = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(mt, seed);

	sigma = 1.0;
	t = 0.0;
	x = 0.0;
	v = 0.0;

	printf("#seed = %ld\n#t\t\tx\t\tv\t\tn_gauss\n", seed);

	while(t < TMAX)
	{
		n_gauss = gsl_ran_gaussian(mt, sigma);
		x = x + v*dt;
		v = v - GAMMA*v - ALPHA*x*dt + (BETA*sqrt(dt))*n_gauss;
		printf("%lf\t%lf\t%lf\t%lf\n", t, x, v, n_gauss);
		t = t+dt;
	}

	gsl_rng_free(mt);

	return 0;
}
