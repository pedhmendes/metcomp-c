/********************************************************
 *                TURING MODEL STABILITY                *
 *******************************************************/      
/********************************************************
 * 	programa para checar a estabilidade das		*
 * 	constantes do modelo de turing			*
 * 	baseado na teoria presente no livro		*
 * 	do sayama					*
 * 							*
 *      executar o programa inserindo                   *
 *      os parametros na seguinte ordem                 *
 *      a b c d Du Dv  	                                *
 *                                                      *
 *      exemplo de exec                                 *
 *      ./a.out 1. -1. 2. -1.5 0.0001 0.0006	        *
 *      ./a.out a   b  c   d   Du     Dv    	        *
 *******************************************************/

/********************************************************
 *                      INCLUDES                        *
 *******************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

/********************************************************
 *                 GLOBAL VARIABLES                     *
 *******************************************************/
double Du, Dv;
double a, b, c, d;
double h, k;

/********************************************************
 *                  MAIN PROGRAM                        *
 *******************************************************/
int main(int argc, char *argv[ ])
{
	double Det = 0;
	double l_side, r_side;
	a = atof(argv[1]);
	b = atof(argv[2]);
	c = atof(argv[3]);
	d = atof(argv[4]);
	Du = atof(argv[5]);
	Dv = atof(argv[6]);

	Det = (a*d)-(b*c);

	l_side = a*Dv + d*Du;
	r_side = 2*(sqrt(Du*Dv*Det));

	if(l_side > r_side)
	{
		printf("\nO sistema eh nao homogenio\nessa escolha de constantes formara padroes\n\n");
	}
	else
	{
		printf("\nO sistema eh homogenio\nessa escolha de constantes nao formara padroes\n\n");
	}

	return 0;
}

