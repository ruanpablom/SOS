
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
#include "sos.h"

#define FAIL    0
#define min(x,y) ((x<y) ? x : y)
#define COND FUNCTION==9 || FUNCTION==10 || FUNCTION==11 || FUNCTION==12 || FUNCTION==13 || FUNCTION==14 || FUNCTION==15


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

	if (GetParameters(argv) == -1){	//read input file
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

	

	
	
	prepararObjFunc();
	

	for (r=0;r<RUN;r++){	
		//Init population
		initPop();
		bestfo = 0.0;
		best_index = 0;
		//Objective function calculation for each individual
		for (i = 0;i<POP_SIZE;i++){
			for (j = 0;j<DIM;j++)
				U[j] = pop[i][j];
	
			fo[i] = objfunc(U, 0);
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
				mutualism_phase(i);
				num_fit_eval+=2;
				commensalism_phase(i);
				num_fit_eval++;
				parasitism_phase(i);
				num_fit_eval++;
				for(j=0;j<POP_SIZE;j++){
					if(fo[j]<=bestfo ){
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
			objfunc(best, 1);
			//values of constraints
			switch(FUNCTION){
				case 9: //Cantilever Beam
					fprintf(file,"g1=%g ",y[0]);
					if(y[0]>1) fprintf(file, "Fail\n");
					else fprintf(file, "Ok\n");
					printf("g1=%g ",y[0]);
					if(y[0]>1)printf("Fail\n");
					else printf("Ok\n");
					break;
				case 10: //I-Beam vertical deflection 
					fprintf(file, "g1=%g ",y[0]);
					if(y[0]>300) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file, "g2=%g ",y[1]);
					if(y[1]>56) fprintf(file, "Fail\n");
					else fprintf(file, "Ok\n");
					printf("g1=%g ",y[0]);
					if(y[0]>300) printf("Fail ");
					else printf("Ok ");
					printf("g2=%g ",y[1]);
					if(y[1]>56) printf("Fail\n");
					else printf("Ok\n");
					break;
				case 11: //Welded Beam 
					fprintf(file,"g1=%g ",y[0]);fprintf(file,"g1=%g ",y[0]);
					if(y[0]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g2=%g ",y[1]);
					if(y[1]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g3=%g ",y[2]);
					if(y[2]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g4=%g ",y[3]);
					if(y[3]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g5=%g ",y[4]);
					if(y[4]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g6=%g ",y[5]);
					if(y[5]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g7=%g ",y[6]);
					if(y[6]>0) fprintf(file, "Fail\n");
					else fprintf(file, "Ok\n");

					printf("g1=%g ",y[0]);
					if(y[0]>0) printf("Fail ");
					else printf("Ok ");
					printf("g2=%g ",y[1]);
					if(y[1]>0) printf("Fail ");
					else printf("Ok ");
					printf("g3=%g ",y[2]);
					if(y[2]>0) printf("Fail ");
					else printf("Ok ");
					printf("g4=%g ",y[3]);
					if(y[3]>0) printf("Fail ");
					else printf("Ok ");
					printf("g5=%g ",y[4]);
					if(y[4]>0) printf("Fail ");
					else printf("Ok ");
					printf("g6=%g ",y[5]);
					if(y[5]>0) printf("Fail ");
					else printf("Ok ");
					printf("g7=%g ",y[6]);
					if(y[6]>0) printf("Fail\n");
					else printf("Ok\n");
					break;
				case 12: //Pressure Vessel 
					fprintf(file,"g1=%g ",y[0]);
					if(y[0]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g2=%g ",y[1]);
					if(y[1]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g3=%g ",y[2]);
					if(y[2]>0) fprintf(file, "Fail ");
					else fprintf(file, "Ok ");
					fprintf(file,"g4=%g ",y[3]);
					if(y[3]>0) fprintf(file, "Fail\n");
					else fprintf(file, "Ok\n");

					printf("g1=%g ",y[0]);
					if(y[0]>0) printf("Fail ");
					else printf("Ok ");
					printf("g2=%g ",y[1]);
					if(y[1]>0) printf("Fail ");
					else printf("Ok ");
					printf("g3=%g ",y[2]);
					if(y[2]>0) printf("Fail ");
					else printf("Ok ");
					printf("g4=%g ",y[3]);
					if(y[3]>0) printf("Fail\n");
					else printf("Ok\n");
					break;
				case 13: //Tension/compression string
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[3],y[3]);
					printf("g1=%g g2=%g g3=%g g4=%g\n",y[0],y[1],y[3],y[3]);
					break;
				case 14: //Speed Reducer(Gear Train)
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10]);
					printf("g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10]);
					break;
				case 15: //10-Bar-Truss
					printf("Stress Violation: ");
					for(i=0;i<10;i++){
						printf("g(%i)=%g ",i+1, y[i]);
					}
					printf("\n");
					printf("DispX: ");
					for(i=10;i<16;i++){
						printf("g(%i)=%g ",i+1, y[i]);
					}
					printf("\n");
					printf("DispY: ");
					for(i=16;i<22;i++){
						printf("g(%i)=%g ",i+1, y[i]);
					}
					printf("\n");
					break;
				/*case 16: //Heat Exchanger Design 
					fprintf(file,"g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7]);
					printf("g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g\n",y[0],y[1],y[3],y[3],y[4],y[5],y[6],y[7]);
					break;*/
			}
			//
		}
		printf("N_fit_eval:");
		fprintf(file,"N_fit_eval:");
		printf("%i \n\n",num_fit_eval);
		fprintf(file, "%i \n\n",num_fit_eval);
		fclose(file);
		//
		if(constr(best)==0)var[r]=bestfo;
		else var[r]=2147483646;

		strf[((strchr(strf,'.'))-strf)]='\0';
		plot(shellComands,strf);//plotting for each run
		
	}//end FOR RUN
	bestfo=var[0];
	best_index=0;
	for (i = 0; i < POP_SIZE; i++){
		if(var[i]<=bestfo && var[i]!=2147483646){
			bestfo=var[i];
			best_index=i;
		}
	}
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
	
	int nfeasible = AvgStdDev(&avg,&stdDev,var);
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
	printf("Feasible: ");
	fprintf(file,"Feasible: ");
	printf("%i\n",nfeasible);
	fprintf(file, "%i\n",nfeasible);
	printf("====================\n");
	fprintf(file,"====================\n");
	fclose(file);
	freeArrays(POP_SIZE, pop, fo, best, ub, lb, y, c);
	free(U);
	
	return 0;
}










