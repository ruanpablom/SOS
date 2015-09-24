#include <stdio.h>
#include <math.h>
#include "mersenne.h"
#include <time.h>
#include <stdlib.h>

//#define RAND_MAX 3

double randon( double inferior, double superior);



int main(){
	int i;
	double sol[10];
	            
	sol[0]=6.01642;
	sol[1]=5.30965;
	sol[2]=4.49434;
	sol[3]=3.50064;
	sol[4]=2.15261;
	sol[5]=197.717;
	sol[6]=197.717;
	sol[7]=197.717;
	sol[8]=197.717;
	sol[9]=197.717;
	printf("%g\n",(0.0624*(sol[0]+sol[1]+sol[2]+sol[3]+sol[4])));
	/*for(i=0;i<10000;i++)
		printf("%g\n",randon(1,4));
	*/
	return 0;
}

double randon( double inferior, double superior){
	double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
	return aux;
}