#include <stdio.h>
#include <math.h>
#include "mersenne.h"
#include <time.h>
#include <stdlib.h>

//#define RAND_MAX 3

double randon( double inferior, double superior);



int main(){
	int i;
	double sol[4];
	    
	sol[0]=0.781374;
	sol[1]=0.386296;
	sol[2]=40.4856;
	sol[3]=197.717;
	sol[4]=197.717;
	sol[5]=197.717;
	sol[6]=197.717;
	sol[7]=197.717;
	sol[8]=197.717;
	sol[9]=197.717;
	//printf("%g\n",((0.6224*sol[0]*sol[2]*sol[3])+(1.7781*sol[1]*pow(sol[2],2))+(3.1661*pow(sol[0],2)*sol[3])+(19.84*pow(sol[0],2)*sol[2])));
	for(i=0;i<10000;i++)
		printf("%g\n",randon(1,4));
	return 0;
}

double randon( double inferior, double superior){
	double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
	return aux;
}