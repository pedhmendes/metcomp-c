//////////////////////////////
//////MOLECULAR DYNAMICS//////
////////PEDRO MENDES//////////
//////////////////////////////

// compile   -> gcc molecular_dynamics1.c -lm
// run	     -> ./a.out | gnuplot

// options:
// rand seed -> gcc -DRAND md5.c -lm
// solid wall-> gcc -DSOLID md5.c -lm

/*******INCLUDES*******/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/*******DEFINES*******/
#define LX 10
#define LY 10
#define NX 5
#define NY 5
#define NP (NX*NY)
#define dt 0.01
#define RC 4
#define RC2 (RC*RC)
#define tmax 50

/*******GLOBAL VARIABLES*******/
double x[NP], y[NP];
double xo[NP], yo[NP];
double vx[NP], vy[NP];
double fx[NP], fy[NP];
double fx_aux[NP], fy_aux[NP];
double Up, Ek, T0=1.0;

/*******FUNCTIONS*******/
void initialize_system(void);
void compute_forces(void);
void velocity_verlet(void);
void gnu_plot(void);
void simulate(void);

/*******MAIN PROGRAM*******/
int main()
{
	simulate();
	return 0;
}

void initialize_system(void)
{
	int i;
	double vcx, vcy, v2, fs;
	double dx, dy;

/*******INITIAL VELOCITIES*******/
	vcx = 0.0; 	vcy = 0.0;
	v2 = 0.0;

	for(i=0; i<NP; i++)
	{
		vx[i] = 1.*rand()/RAND_MAX;
		vy[i] = 1.*rand()/RAND_MAX;

		vcx += vx[i];
		vcy += vy[i];

		v2 += vx[i]*vx[i] + vy[i]*vy[i];
	}

	vcx /= NP; 	vcy /= NP;
	v2 /= NP;

	fs = sqrt(2.*T0/v2);

	for(i=0; i<NP; i++)
	{
		vx[i] = (vx[i]-vcx)*fs;
		vy[i] = (vy[i]-vcy)*fs;
	}

/*******INITIAL POSITIONS*******/
	dx = (double)LX/NX;
	dy = (double)LY/NY;

	for(i=0; i<NP; i++)
	{
		x[i] = (i%NX)*dx;
		y[i] = (i/NY)*dy;

		xo[i] = x[i] - vx[i]*dt;
		yo[i] = y[i] - vy[i]*dt;
	}

	return;
}
	
void compute_forces(void)
{
	int i, j;
	double x1, x2, dx;
	double y1, y2, dy;
	double r2, r2i, r6, ff;
	
	Up = 0.0;

/*******CLEAN ARRAYS*******/
	for(i=0; i<NP; i++)
	{
		fx[i] = 0.0;
		fy[i] = 0.0;
	}

/*******FOR ALL PARTICLES*******/
	for(i=0; i<(NP-1); i++)
	{
		x1 = x[i];
		y1 = y[i];

		for(j=i+1; j<NP; j++)
		{
			x2 = x[j];
			dx = x1-x2;
			dx = fmod(dx,LX);
			dx = dx - rint(dx/LX)*LX;

			y2 = y[j];
			dy = y1-y2;
			dy = fmod(dy,LY);
			dy = dy - rint(dy/LY)*LY;

			r2 = dx*dx + dy*dy;

			if(r2 < RC2)
			{
				r2i = 1/r2;
				r6 = r2i*r2i*r2i;
				ff = 48*r2i*r6*(r6-0.5);

				fx[i] += ff*dx; 	fy[i] += ff*dy;
				fx[j] -= ff*dx;		fy[j] -= ff*dy;
		
				Up += 4*r6*(r6-1);
			}
		}
	}

	return;
}
void velocity_verlet(void)
{
/*******INTEGRATE EQUATIONS*******/
	int i;

	for(i=0; i<NP; i++)
	{
		x[i] = x[i] + vx[i]*dt + fx[i]*dt*dt*0.5;
		y[i] = y[i] + vy[i]*dt + fy[i]*dt*dt*0.5;

		fx_aux[i] = fx[i];
		fy_aux[i] = fy[i];
	}

	compute_forces();

	for(i=0; i<NP; i++)
	{
		vx[i] = vx[i] + 0.5*(fx_aux[i]+fx[i])*dt;
		vy[i] = vy[i] + 0.5*(fy_aux[i]+fy[i])*dt;
	}

/*******SOLID WALLS*******/
#ifdef SOLID
	for(i=0; i<NP; i++)
	{
		if(x[i] > LX)
		{
			x[i] = 2*LX-x[i];
			vx[i] = -vx[i];
		}
		if(x[i] < 0.)
		{
			x[i] = -x[i];
			vx[i] = -vx[i];
		}
		if(y[i] > LY)
		{
			y[i] = 2*LY-y[i];
			vy[i] = -vy[i];
		}
		if(y[i] < 0.)
		{
			y[i] = -y[i];
			vy[i] = -vy[i];
		}
	}
#endif

	return;
}

void gnu_plot(void)
{
/*******GNUPLOT VISUALIZATION*******/
	int i;

	printf("set xrange [0:%d]\nset yrange [0:%d]\n", LX,LY);
	printf("set size square\n");
	printf("set key off\n");
	printf("unset xtics\nunset ytics\n");
	printf("unset colorbox\n");
	printf("set style circle radius (0.5*(2**(1./6.)))\n");
	printf("set style fill solid 0.6\n");
	printf("plot \"-\" u 2:3 w circles lc rgb 'blue'\n");
	
	for(i=0; i<NP; i++)
	{
		printf("%d\t%lf\t%lf\n", i, fmod(fmod(x[i],LX)+LX,LX), fmod(fmod(y[i],LY)+LY,LY));
	}

	printf("e\n\n");

	return;
}

void simulate(void)
{
/*******SIMULATION*******/
	int i;
	double t=0;
	char F_1[100];
	FILE *f_1;

#ifdef RAND
	srand(time(NULL));
#endif

	sprintf(F_1, "data_NP%d_T%.2lf.dat", NP, T0); 
	f_1 = fopen(F_1, "w");

	initialize_system();
	compute_forces();

/*******REFERENCE POSITIONS*******/
	for(i=0; i<NP; i++)
	{
		xo[i] = x[i];
		yo[i] = y[i];
	}

	fprintf(f_1, "#t\tUg\tEk\tEM");

/*******SAMPLING*******/
	while(t <= tmax)
	{
		velocity_verlet();

		printf("set title 't = %.3lf'\n", t);
		gnu_plot();

/*******KINETIC ENERGY*******/
		Ek = 0.0;
	
		for(i=0; i<NP; i++)
		{
			Ek += ((vx[i]*vx[i])+(vy[i]*vy[i]))*0.5;
		}

		fprintf(f_1, "%lf\t%lf\t%lf\t%lf\n", t, Up, Ek, (Up+Ek));

		t = t + dt;
	}

	printf("pause 2\n");

	return;
}
