
/* 
	Description: The ANSI C code of the SOS Algorithm
	Programmer:  Ruan Pablo Medeiros
	E-mail:      pm.ruan@gmail.com
	Date:	     04/09/2014
	Lisence:     Free
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "mersenne.h"
#include "sosfunctions.h"

#define FAIL    0
#define min(x,y) ((x<y) ? x : y)
#define COND FUNCTION!=10 && FUNCTION!=11 && FUNCTION!=12 && FUNCTION!=13 && FUNCTION!=14

/* Control Parameters of the search algorithm*/
int POP_SIZE;  /* The number of candidate solutions*/
int MAX_ITER; /*The number of iterations*/
int FUNCTION;
//Problem definitions
int DIM;    //number of problem variables
double *lb; //lower bound of the variables
double *ub; //upper bound of the variables
int RUN;  /*Algorithm can run many times in order to see its robustness*/
double bestfoRUN;

//Global variables
	double **pop; //[POP_SIZE][DIM];  //population of candidate solutions.
	double *fo; //[POP_SIZE];      //objective function value.
	double *best; //[DIM];           //best solution found
	double bestfo;         //best fo value
	int best_index;     	  //index for the best solution
	double *y;

void AllocArrays();

/*Main program of the search algorithm*/
int main(int argc, char **argv){
	
	int i, j, k, r;
	double *U; //[DIM]; // vetor para cálculo da func obj.
	double avg;
	double stdDev;
	double *var;
	int num_fit_eval=0;
	//int max_fit_eval;
	int num_iter=0;
	double mediafo=0;
	char str[]="dadosplot//dadosplot";
	char strf[100];
	double *mediaBfo;
	double *mediaM;
	FILE *file;
	FILE *shellComands;
	

	//refresh the folder of plotting
	shellComands = popen ("rm dadosplot// -R", "w");
	pclose(shellComands);
	shellComands = popen ("mkdir dadosplot", "w");
	pclose(shellComands);
	//

	srand(time(NULL));
	MT_seed();

	if (GetParameters(argv, &RUN, &MAX_ITER, &POP_SIZE, &DIM, &FUNCTION) == -1){	//read input file
		return 0;
	}
	showParameters(FUNCTION, RUN, MAX_ITER, POP_SIZE, DIM);
	
	//Alloc arrays
	mediaBfo=(double*)malloc(MAX_ITER*sizeof(double));
	mediaM=(double*)malloc(MAX_ITER*sizeof(double));
	for(i=0;i<MAX_ITER;i++){
		mediaBfo[i]=0;
		mediaM[i]=0;
	}
	var=(double*)malloc(RUN*sizeof(double));
	


	AllocArrays();	

	U = (double*)malloc(DIM * sizeof(double));
	//

	

	
	
	prepararObjFunc(FUNCTION, lb, ub);
	

	for (r=0;r<RUN;r++){	
		//Init population
		initPop(FUNCTION, POP_SIZE, DIM, best, pop, lb, ub, fo);
		bestfo = 0.0;
		best_index = 0;
		//Objective function calculation for each individual
		for (i = 0;i<POP_SIZE;i++){
			for (j = 0;j<DIM;j++)
				U[j] = pop[i][j];
	
			fo[i] = objfunc(U, y, FUNCTION, DIM);
		}
		//Best current solution identification.
		bestfo = fo[0];
		best_index = 0;
		for(i=0;i<POP_SIZE;i++){
			if (fo[i]<=bestfo) {
	       		bestfo=fo[i];
	       		for(j=0;j<DIM;j++){
	        	   best[j]=pop[i][j];
				}
			best_index = i;
	        }
		}

		num_fit_eval=0;
		//max_fit_eval;
		num_iter=0;
		mediafo=0;
		strcpy(str,"dadosplot//dadosplot");
		converteDecChar(strf,r);
		strcat(strf,".txt");
		strcat(str,strf);
		if((file = fopen(str,"a")) == NULL){
    	    	printf("Erro ao abrir arquivo!!!\n\n");
    	    	exit(1);
    	  	}
		fprintf(file,"%s %14s %15s\n","#ITER","#BEST_FO","#MEDIA_FO");
		fclose(file);
		while(num_iter<MAX_ITER){
			num_iter++;
			for(i=0;i<POP_SIZE;i++){
				mutualism_phase(y, i, DIM, pop, best, ub, lb, fo, FUNCTION, POP_SIZE);
				num_fit_eval+=2;
				commensalism_phase(y, i, DIM, pop, best, ub, lb, fo, FUNCTION, POP_SIZE);
				num_fit_eval++;
				parasitism_phase(y, i, DIM, pop, ub, lb, fo, FUNCTION, POP_SIZE);
				num_fit_eval++;
				for(j=0;j<POP_SIZE;j++){
					if(fo[j]<=bestfo){
						bestfo=fo[j];
						best_index=j;
					}
				}
				for(j=0;j<DIM;j++){
					best[j]=pop[best_index][j];
				}	
			}
			for(i=0;i<POP_SIZE;i++){
				mediafo+=fo[i];
			}
			mediafo=mediafo/POP_SIZE;
			if((file = fopen(str,"a")) == NULL){
       			printf("Erro ao abrir arquivo!!!\n\n");
        		exit(1);
      		}
			fprintf(file,"%i%18g%15g\n",num_iter,bestfo,mediafo);
			fclose(file);
			mediaM[num_iter]+=mediafo;//sum of all mediafo in the num_iter position
			mediaBfo[num_iter]+=bestfo;//sum of all bestfo	in the num_iter position
		}
	
		//Loop de Iterações.
	
		if((file = fopen("dadosplot//exec.txt","a")) == NULL){
       		printf("Erro ao abrir arquivo!!!\n\n");
       		exit(1);
      	}
		printf("RUN: %d\n",r);
		fprintf(file,"RUN: %d\n",r);
		printf("Best solution: ");
		fprintf(file,"Best solution: ");
		for (k=0; k<DIM;k++){//variables
			printf("%g ",best[k]);
			fprintf(file,"%g ",best[k]);
		}
		printf(" Fo:");
		fprintf(file," Fo:");
		printf("%g \n",bestfo);
		fprintf(file,"%g \n",bestfoRUN);
		if(r==0)bestfoRUN=bestfo;
		bestfoRUN=min(bestfo, bestfoRUN);
		printf("MIN: %g\n",bestfoRUN);
		printf("bestfo: %g\n", bestfo);
		printf("bestfoRUN: %g\n", bestfoRUN);
	
		if(COND){
			objfunc(best, y, FUNCTION, DIM);
			//values of constraints
			switch(FUNCTION){
				case 9: //Cantilever Beam
					fprintf(file,"g1=%g\n",y[0]);
					printf("g1=%g\n",y[0]);
					break;
				case 10: //I-Beam vertical deflection 
					fprintf(file,"g1=%g g2=%g\n",y[0],y[1]);
					printf("g1=%g g2=%g\n",y[0],y[1]);
					break;
				case 11: //Welded Beam 
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6]);
					printf("g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6]);
					break;
				case 12: //Pressure Vessel 
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[2],y[3]);
					printf("g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[2],y[3]);
					break;
				case 13: //Tension/compression string
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[3],y[3]);
					printf("g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[3],y[3]);
					break;
				case 14: //Speed Reducer(Gear Train)
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10]);
					printf("g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10]);
					break;
			}
			//
		}
		printf("N_fit_eval:");
		fprintf(file,"N_fit_eval:");
		printf("%i \n\n",num_fit_eval);
		fprintf(file, "%i \n\n",num_fit_eval);
		fclose(file);
		//
		var[r]=bestfo;
		strf[((strchr(strf,'.'))-strf)]='\0';
		plot(shellComands,strf);//plotting for each run
		
	}//end FOR RUN
	
	if((file = fopen("dadosplot//dadosplotFinal.txt","a")) == NULL){
        printf("Erro ao abrir arquivo!!!\n\n");
        exit(1);
    }
	fprintf(file,"%s %14s %15s\n","#ITER","#MEDIA_BFO","#MEDIA_M");
	for(i=1;i<=MAX_ITER;i++){
		fprintf(file,"%i%18g%15g\n",i,(mediaBfo[i]/r),(mediaM[i]/r));
	}
	fclose(file);
	strcpy(strf,"Final");
	plot(shellComands,strf);//final plotting
	
	if((file = fopen("dadosplot//exec.txt","a")) == NULL){
       	printf("Erro ao abrir arquivo!!!\n\n");
       	exit(1);
    }
	AvgStdDev(&avg,&stdDev,var,RUN);
	printf("====================\n");
	fprintf(file,"====================\n");
	printf("Best Fo: ");
	fprintf(file,"Best Fo: ");
	printf("%g\n",bestfoRUN);
	fprintf(file,"%g\n",bestfoRUN);
	printf("Avg: ");
	fprintf(file,"Avg: ");
	printf("%g\n",avg);
	fprintf(file,"%g\n",avg);
	printf("StdDev: ");
	fprintf(file,"StdDev: ");
	printf("%g\n",stdDev);
	fprintf(file, "%g\n",stdDev);
	printf("====================\n");
	fprintf(file,"====================\n");
	fclose(file);
	freeArrays(POP_SIZE, pop, fo, best, ub, lb);
	free(U);
	
	return 0;
}

void AllocArrays(){
	int j;
	int constn=0;

	pop = (double**)malloc (POP_SIZE*sizeof(double*));
        for (j = 0;j < POP_SIZE; j++)
		pop[j] = (double*)malloc (DIM * sizeof(double));

	fo = (double*)malloc (POP_SIZE*sizeof(double));

	best = (double*)malloc (DIM*sizeof(double));

	ub=(double*)malloc (DIM*sizeof(double));
	lb=(double*)malloc (DIM*sizeof(double));

	switch(FUNCTION){
		case 9:
			y=(double*)malloc(sizeof(double));
			break;
		case 10:
			constn=2;
			y=(double*)malloc(constn*sizeof(double));
			break;
		case 11:
			constn=7;
			y=(double*)malloc(constn*sizeof(double));
			break;
		case 12: //Pressure Vessel
			constn=4;
			y=(double*)malloc(constn*sizeof(double));
			
			break;
		case 13: //Tension/compression string
			constn=4;
			y=(double*)malloc(constn*sizeof(double));
			break;
		case 14: // Speed Reducer(Gear Train)
			constn=11;
			break;
	}

}








