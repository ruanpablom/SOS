
/* 
Description: The ANSI C code of the SOS Algorithm
Programmer:  Ruan Pablo Medeiros
E-mail:      pm.ruan@gmail.com
Date:	     04/09/2014
Lisence:     Free
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
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
    int i, k, r;
    double avg;
    double stdDev;
    double *var;
    double mediafo=0;
    char str[]="dadosplot//dadosplot";
    char strf[100];
    double *mediaBfo;
    double *mediaM;
    FILE *file;
    FILE *shellComands;

/*    //refresh the folder of plotting
    shellComands = popen ("rm dadosplot// -R", "w");
    pclose(shellComands);
    shellComands = popen ("mkdir dadosplot", "w");
    pclose(shellComands);
    //
*/
    srand(time(NULL));

    if (GetParameters(argv) == -1){	//read input file
        return 0;
    }
    showParameters(FUNCTION, RUN, MAX_ITER, POP_SIZE, DIM);

    switch(FUNCTION){
        case 9:
            rest=1;
            break;
        case 10:
            rest=2;
            break;
        case 11:
            rest=7;
            break;
        case 12: //Pressure Vessel
            rest=4;
            break;
        case 13: //Tension/compression string
            rest=4;
            break;
        case 14: // Speed Reducer(Gear Train)
            rest=11;
            break;
        case 15: // 10-Bar Truss
            rest=21;
            break;
    }

    //Alloc arrays
    mediaBfo=(double*)malloc(MAX_ITER*sizeof(double));
    mediaM=(double*)malloc(MAX_ITER*sizeof(double));
    for(i=0;i<MAX_ITER;i++){
        mediaBfo[i]=0;
        mediaM[i]=0;
    }

    var=(double*)malloc(RUN*sizeof(double));

    AllocArrays();	

    prepararObjFunc();

    for (r=0;r<RUN;r++){	
        //Init population
        initPop();
        bestfo = 10000000;

        for(i=0;i<POP_SIZE;i++){
            if(fo[i]<=bestfo){
                best_index=i;
                bestfo=fo[i];
            }
        }
        for(i=0;i<DIM;i++)best[i]=pop[best_index][i];

        num_fit_eval=0;
        num_iter=0;
        mediafo=0;
        /* strcpy(str,"dadosplot//dadosplot");
           converteDecChar(strf,r);
           strcat(strf,".txt");
           strcat(str,strf);
           if((file = fopen(str,"a")) == NULL){
           printf("Erro ao abrir arquivo!!!\n\n");
           exit(1);
           }
           fprintf(file,"%s %14s %15s\n","#ITER","#BEST_FO","#MEDIA_FO");
           fclose(file);*/
        while(num_iter<MAX_ITER){
            num_iter++;
            sos_iter();
            for(i=0;i<POP_SIZE;i++){
                mediafo+=fo[i];
            }
            mediafo=mediafo/POP_SIZE;
            /*if((file = fopen(str,"a")) == NULL){
              printf("Erro ao abrir arquivo!!!\n\n");
              exit(1);
              }
              fprintf(file,"%i%18g%15g\n",num_iter,bestfo,mediafo);
              fclose(file);*/
            mediaM[num_iter]+=mediafo;//sum of all mediafo in the num_iter position
            mediaBfo[num_iter]+=bestfo;//sum of all bestfo	in the num_iter position
        }

        //Loop de Iterações.
        /*
           if((file = fopen("dadosplot//exec.txt","a")) == NULL){
           printf("Erro ao abrir arquivo!!!\n\n");
           exit(1);
           }*/
        /*printf("RUN: %d\n",r);
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

          showConst(var,r,file);



          printf("N_fit_eval:");
          fprintf(file,"N_fit_eval:");
          printf("%i \n\n",num_fit_eval);
          fprintf(file, "%i \n\n",num_fit_eval);
          fclose(file);*/
        //
        /*
           strf[((strchr(strf,'.'))-strf)]='\0';
           plot(shellComands,strf);//plotting for each run*/

    }//end FOR RUN
    /* bestfo=var[0];
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
    //	freeArrays();
    */
    return 0;
}
