///////////////////////////
////////PEDRO MENDES///////
///////////////////////////

// compile	-> gcc -Wall ho_baoab_gsl.c -lgsl -lgslcblas -lm
// run		-> ./a.out > serie.dsf
	
#include<stdio.h>
#include<time.h>
#include<math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>

#define 	TMAX		100 
#define 	DT		0.001
#define		GAMMA		0.5
#define 	M		1
#define		KB		1
#define		TEMP		2

int main()
{
	double t, x, p;
	double f, k;
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
	k = 1.0;
	f = -k*x;

	printf("#seed = %ld\n#t\t\tx\t\tp\t\tf\n", seed);

	while(t < TMAX)
	{
		n_gauss = gsl_ran_gaussian(mt, sigma);

		p = p + DT*0.5*f;
		x = x + DT*0.5*(p/M);
		p = exp(-GAMMA*DT)*p + sqrt(1-exp(-2*GAMMA*DT))*sqrt(KB*M*TEMP)*n_gauss;
		x = x + DT*0.5*(p/M);
		f = -k*x;
		p = p + DT*0.5*f;

		printf("%lf\t%lf\t%lf\t%lf\n", t, x, p, f);
		t = t + DT;
	}

	gsl_rng_free(mt);

	return 0;
}
