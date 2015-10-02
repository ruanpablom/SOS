
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include "sos.h"
#include "mersenne.h"

#define FAIL 0
#define PIII 3.14159265359
#define COND FUNCTION==9 || FUNCTION==10 || FUNCTION==11 || FUNCTION==12 || FUNCTION==13 || FUNCTION==14 || FUNCTION==15

pthread_mutex_t data_mutex = PTHREAD_MUTEX_INITIALIZER;

double randon( double inferior, double superior){
	double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
	return aux;
}

void freeArrays(){
	free(fo);
	free(best);
	free(lb);
	free(ub);
	free(pop[0]);	
	free(pop);
	free_slice();
}

int GetParameters(char **argv){
	FILE *file = fopen( argv[1], "r" );
    
	if (file == 0){
    	printf( "Could not open ini file! Usage ./<exec> <file.in>\n" );
		return -1;
	}else{
		ffscanf("RUN", file, "%d", &RUN);
		ffscanf("MAX_ITER", file, "%d", &MAX_ITER);
		ffscanf("POP_SIZE", file, "%d", &POP_SIZE);
		ffscanf("DIM", file, "%d", &DIM); 
		ffscanf("FUNCTION",file, "%d", &FUNCTION);
		ffscanf("CORES",file, "%d", &CORES);
	return 1;
    }
    fclose(file);
}

void showParameters(){
	printf("***PARAMETERS***\n");
	printf("RUNS = %d\n", RUN);
	printf("MAX_ITER = %d\n", MAX_ITER);
	printf("POP_SIZE = %d\n", POP_SIZE);
	printf("DIM = %d\n", DIM);
	printf("CORES = %d\n", CORES);
	switch (FUNCTION){
		case 0:
			printf("FUNCTION = %s\n","Rastrigin");
			break;
		case 1:
			printf("FUNCTION = %s\n","Schaffer");
			break;
		case 2:
			printf("FUNCTION = %s\n","Griewank");
			break;
		case 3:
			printf("FUNCTION = %s\n","Ackley");
			break;
		case 4:
			printf("FUNCTION = %s\n","Rosenbrock");
			break;
		case 5:
			printf("FUNCTION = %s\n","Sphere");
			break;
		case 6:
			printf("FUNCTION = %s\n","Michaelewicz10");
			break;
		case 7:
			printf("FUNCTION = %s\n","Booth");
			break;
		case 8:
			printf("FUNCTION = %s\n","Quartic");
			break;
		case 9:
			printf("FUNCTION = %s\n","Cantilever beam");
			break;
		case 10:
			printf("FUNCTION = %s\n","I-Beam vertical deflection");
			break;
		case 11:
			printf("FUNCTION = %s\n","Welded Beam");
			break;
		case 12:
			printf("FUNCTION = %s\n","Pressure Vessel");
			break;
		case 13:
			printf("FUNCTION = %s\n","Tension/compression string");
			break;
		case 14:
			printf("FUNCTION = %s\n","Speed Reducer(Gear Train)");
			break;
		case 15:
			printf("FUNCTION = %s\n","10-bar truss");
			break;
	}
	printf("****************\n");
}

void prepararObjFunc(){
    switch (FUNCTION){
    	case 0: //Rastrigin
        	lb[0] = -5.12;
        	ub[0] = 5.12;
        	break;
    	case 1: //Schaffer
        	lb[0] = -100.00;
        	ub[0] = 100.00;
        	break;
    	case 2: //Griewank
       		lb[0] = -600.00;
        	ub[0] = 600.00;
        	break;
    	case 3: //Ackley
        	lb[0] = -32.00;
        	ub[0] = 32.00;
        	break;
    	case 4: //Rosenbrock
        	lb[0] = -30.00;
        	ub[0] = 30.00;
        	break;
    	case 5: //Sphere
			lb[0] = -100.00;
			ub[0] = 100.00;
			break;
    	case 6: //Michaelbewicz10
			lb[0] = 0.00;
			ub[0] = M_PI;
			break;
    	case 7: //Booth
			lb[0] = -10.00;
			ub[0] = 10.00;
			break;
    	case 8: //Quartic
			lb[0] = -1.28;
			ub[0] = 1.28;
			break;
		case 9: //Cantilever beam
			lb[0] = 0.01;
			ub[0] = 100.0;
			break;
		case 10: //10: I-Beam vertical deflection 
			//sol[0]=b
    		//sol[1]=h
    		//sol[2]=tw
    		//sol[3]=tf
			//b domain
			lb[0] = 10.0;
			ub[0] = 50.0;
			//
			//h domain
			lb[1] = 10.0;
			ub[1] = 80.0;
			//
			//tw domain
			lb[2] = 0.9;
			ub[2] = 5.0;
			//
			//tf domain
			lb[3] = 0.9;
			ub[3] = 5.0;
			//
			break;
		case 11: // Welded Beam 
			lb[0] = 0.1;
			ub[0] = 2.0;
			lb[1] = 0.1;
			ub[1] = 10.0;
			lb[2] = 0.1;
			ub[2] = 10.0;
			lb[3] = 0.1;
			ub[3] = 2.0;
			//
			break;
		case 12: // Pressure Vessel
			lb[0] = 0.0625;
			ub[0] = 99.0*0.0625;
			lb[1] = 0.0625;
			ub[1] = 99.0*0.0625;
			lb[2] = 10.0;
			ub[2] = 200.0;
			lb[3] = 10.0;
			ub[3] = 200.0;
			//
			break;
		case 13: // Tension/compression string
			lb[0] = 0.05;
			ub[0] = 2.0;
			lb[1] = 0.25;
			ub[1] = 1.3;
			lb[2] = 2.0;
			ub[2] = 15.0;
			//
			break;
		case 14: // Speed Reducer(Gear Train)
			lb[0] = 2.6;
			ub[0] = 3.6;
			lb[1] = 0.7;
			ub[1] = 0.8;
			lb[2] = 17;
			ub[2] = 28;
			lb[3] = 7.3;
			ub[3] = 8.3;
			lb[4] = 7.3;
			ub[4] = 8.3;
			lb[5] = 2.9;
			ub[5] = 3.9;
			lb[6] = 5.0;
			ub[6] = 5.5;
			//
			break;
		case 15:
			lb[0] = 0.10;
			ub[0] = 35.00;
			break;
    	default:
    	    printf("Info: Invalid function\n") ;
    	    exit(0);
    	}
}	

double constr(double *sol, int cond){//calculate penalization
	double y[rest];
	double pen=0;
	double r=0,m_1=0,j_1=0,s=0,d=0,pc=0,t1=0,t2=0,t=0,p=0,l=0,e=0,g_1=0,t_max=0,s_max=0,d_max=0;

	double stress[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	double displX[6] = {0., 0., 0., 0., 0., 0.};
	double displY[6] = {0., 0., 0., 0., 0., 0.};

	double P  = 10000.;

	double E  = 10000.;

	double *area;

	int num_node = 0, 
    num_elem = 0;
	int kode[6] = {0, 0, 0, 0, 0, 0};
	double coordx[6] = {0., 0., 0., 0., 0., 0.},
    coordy[6] = {0., 0., 0., 0., 0., 0.},
    fX[6] = {0., 0., 0., 0., 0., 0.},
    fY[6] = {0., 0., 0., 0., 0., 0.};
	int lm[4] = {0,0,0,0}, 
    neq = 0;
	int member[10][2];
	double ymod[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	double a[12][12], 
    b[12] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.}, 
    stiffmatrix[4][4];

	int n = 0, nn = 0;

	double dx = 0.,dy = 0.,xl = 0.;
	double du = 0.,dv = 0.,dl = 0.;

	int nl = 0;
	int n1 = 0,m1 = 0;
	int neq1 = 0;

	int k = 0, i=0, j=0, g=0, m=0;
	double cosa = 0.,sina = 0.,comm = 0.;
	switch(FUNCTION){
		case 9:
			y[0]=((61/pow(sol[0],3))+(37/pow(sol[1],3))+(19/pow(sol[2],3))+(7/pow(sol[3],3))+(1/pow(sol[4],3)));
			if(y[0]>1)pen+=pow(y[0],2);
			pen*=1000;
			break;
		case 10:
			//sol[0]=b
    		//sol[1]=h
    		//sol[2]=tw
    		//sol[3]=tf
			y[0]=2*sol[0]*sol[2]+sol[2]*(sol[1]-(2*sol[3]))+1;
			if(y[0]>300)pen+=y[0]+10000;
			y[1]=(((18*sol[1]*10000)/(sol[2]*pow(sol[1]-(2*sol[3]),3)+2*sol[0]*sol[2]*(4*pow(sol[3],2)+3*sol[1]*(sol[1]-2*sol[3]))))+((15*sol[0]*1000)/(((sol[1]-2*sol[3])*pow(sol[2],3))+2*sol[2]*pow(sol[0],3))));
			if(y[1]>56)pen+=y[1]+10000;
			break;
		case 11://welded beam
			p = 6000; 
			l = 14; 
			e = 30e+6; 
			g_1 = 12e+6;
    		t_max = 13600; 
    		s_max = 30000; 
    		d_max = 0.25;
			j_1=2*(sqrt(2)*sol[0]*sol[1]*(pow(sol[1],2)/12+pow(((sol[0]+sol[2])/2),2)));
			r=sqrt((pow(sol[1],2)/4)+pow((sol[0]+sol[2]/2),2));
			m_1=p*(l+(sol[1]/2));
			t2=(m_1*r)/j_1;
			t1=p/(sqrt(2)*sol[0]*sol[1]);
			t=sqrt(pow(t1,2)+2*t1*t2*(sol[1]/(2*r))+pow(t2,2));
			s=(6*p*l)/(sol[3]*pow(sol[2],2));
			d=(4*p*pow(l,3))/(sol[3]*pow(sol[2],3)*e);
			pc=4.013*e*sqrt((pow(sol[2],2)*pow(sol[3],6))/36)*(1-sol[2]*sqrt(e/(4*g_1))/(2*l))/(l*l);
			y[0]=t-t_max;
			if(y[0]>0)pen+=pow(y[0],2);
			y[1]=s-s_max;
			if(y[1]>0)pen+=pow(y[1],2);
			y[2]=sol[0]-sol[3];
			if(y[2]>0)pen+=pow(y[2],2);
			y[3]=0.10471*pow(sol[0],2)+0.04811*sol[2]*sol[3]*(14+sol[1])-5;
			if(y[3]>0)pen+=pow(y[3],2);
			y[4]=0.125-sol[0];
			if(y[4]>0)pen+=pow(y[4],2);
			y[5]=d-d_max;
			if(y[5]>0)pen+=pow(y[5],2);
			y[6]=6000-pc;
			if(y[6]>0)pen+=pow(y[6],2);
			pen*=1000000;
			break;
		case 12: //Pressure Vessel
			y[0]=-sol[0]+0.0193*sol[2];
			if(y[0]>0){
				pen+=y[0]*y[0];
			}
			y[1]=-sol[1]+0.00954*sol[2];
			if(y[1]>0){
				pen+=y[1]*y[1];
			}
			y[2]=-PIII*pow(sol[2],2)*sol[3]-(4*PIII/3)*pow(sol[2],3)+1296000;
			if(y[2]>0){
				pen+=y[2]*y[2];
			}
			y[3]=sol[3]-240;
			if(y[3]>0){
				pen+=y[3]*y[3];
			}
			pen*=1000000000;
		break;
		case 13: //Tension/compression string
			y[0]=1-((pow(sol[1],3)*sol[2])/(71785*pow(sol[0],4)));
			if(y[0]>0)pow(pen+=y[0],2);
			y[1]=((((4*pow(sol[1],2))-(sol[0]*sol[1]))/(12566*((sol[1]*pow(sol[0],3))-pow(sol[0],4))))+(1/(5108*pow(sol[0],2)))-1);
			if(y[1]>0)pen+=pow(y[1],2);
			y[2]=1-((140.45*sol[0])/(sol[2]*pow(sol[1],2)));
			if(y[2]>0)pen+=pow(y[2],2);
			y[3]=((sol[0]+sol[1])/1.5)-1;
			if(y[3]>0)pen+=pow(y[3],2);
			pen*=1000000;
			break;
		case 14: // Speed Reducer(Gear Train)
			y[0]=27*pow(sol[0],-1)*pow(sol[1],-2)*pow(sol[2],-1)-1;
			if(y[0]>0)pen+=pow(y[0],2);//
			y[1]=397.5*pow(sol[0],-1)*pow(sol[1],-2)*pow(sol[2],-2)-1;
			if(y[1]>0)pen+=pow(y[1],2);
			y[2]=1.93*pow(sol[1],-1)*pow(sol[2],-1)*pow(sol[3], 3)*pow(sol[5],-4)-1;
			if(y[2]>0)pen+=pow(y[2],2);
			y[3]=1.93*pow(sol[1],-1)*pow(sol[2],-1)*pow(sol[4], 3)*pow(sol[6],-4)-1;
			if(y[3]>0)pen+=pow(y[3],2);
			y[4]=(pow((pow(((745*sol[3])/(sol[1]*sol[2])),2)+(16.9*pow(10,6))),0.5)/(0.1*pow(sol[5],3)))-1100;
			if(y[4]>0)pen+=pow(y[4],2);
			y[5]=(pow((pow(((745*sol[4])/(sol[1]*sol[2])),2)+(157.5*pow(10,6))),0.5)/(0.1*pow(sol[6],3)))-850;
			if(y[5]>0)pen+=pow(y[5],2);
			y[6]=sol[1]*sol[2]-40;
			if(y[6]>0)pen+=pow(y[6],2);
			y[7]=(sol[0]/sol[1])-12;
			if(y[7]>0)pen+=pow(y[7],2);
			y[8]=5-(sol[0]/sol[1]);
			if(y[8]>0)pen+=pow(y[8],2);
			y[9]=(1.5*sol[5]+1.9)*pow(sol[3],-1)-1;
			if(y[9]>0)pen+=pow(y[9],2);
			y[10]=(1.1*sol[6]+1.9)*pow(sol[4],-1)-1;
			if(y[10]>0)pen+=pow(y[10],2);
			pen*=1000000000000;
			break;

		case 15: //10-bar-truss
			/*init */
			y[0]=0;
			y[1]=0;
			
			for (i = 0; i < 12; i++){
				for (j = 0; j < 12; j++){
					a[i][j] = 0.;
				}	
			}	

			for (i = 0; i < 4; i++){
				for (j = 0; j < 4; j++){
				       stiffmatrix[i][j] = 0.;
				}
			}

			area = (double*) malloc (DIM * sizeof(double));
			if(area==NULL)exit(0);
			for (i=0; i < DIM; i++){
				area[i] = sol[i];
			}
			
			num_node = 6;
			num_elem = 10;

			kode[4] = 3; /* fixo X e Y*/
			coordx[4] = 0;
			coordy[4] = 0;
			fX[4] = 0;
			fY[4] = 0;
			kode[2] = 0;
			coordx[2] = 360;
			coordy[2] = 0;
			fX[2] = 0;
			fY[2] = 50;
			kode[0] = 0;
			coordx[0] = 720;
			coordy[0] = 0;
			fX[0]=0;
			fY[0] = 50;
			kode[5] = 3;
			coordx[5] = 0;
			coordy[5] = -360;
			fX[5] = 0;
			fY[5] = 0;
			kode[3] = 0;
			coordx[3] = 360;
			coordy[3] = -360;
			fX[3] = 0;
			fY[3] = -150;
			kode[1] = 0;
			coordx[1] = 720;
			coordy[1] = -360;
			fX[1] = 0;
			fY[1] = -150;

			member[0][0] = 5; member[0][1] = 3;
			member[1][0] = 3; member[1][1] = 1;
			member[2][0] = 6; member[2][1] = 4;
			member[3][0] = 4; member[3][1] = 2;
			member[4][0] = 3; member[4][1] = 4;
			member[5][0] = 1; member[5][1] = 2;
			member[6][0] = 5; member[6][1] = 4;
			member[7][0] = 6; member[7][1] = 3;
			member[8][0] = 3; member[8][1] = 2;
			member[9][0] = 1; member[9][1] = 4;

			ymod[0] = E;
			ymod[1] = E;
			ymod[2] = E;
			ymod[3] = E;
			ymod[4] = E;
			ymod[5] = E;
			ymod[6] = E;
			ymod[7] = E;
			ymod[8] = E;
			ymod[9] = E;
			/*end init
			
			//stiff */
			/* Initialize Stiffness Matrix and Load Vector */
						
			neq = 2 * num_node;
			for (m=0; m<neq;m++){
				b[m] = 0.;
				for(k=0;k<neq;k++) a[m][k]=0.;
			}
			/*  >>>>>>> loop over all elements  <<<<<<< */
			for(m=0;m<num_elem;m++){
				/* Determination of Length, Direction */
				i=member[m][0] - 1;
				j=member[m][1] - 1;
				dx = coordx[j]-coordx[i];
				dy = coordy[j]-coordy[i];
				xl = sqrt((dx * dx) + (dy * dy));
				cosa=dx/xl;
				sina=dy/xl;
				comm=area[m]*ymod[m]/xl;
				/* Construction of the 4X4 Element Stiffness Matrix */
				stiffmatrix[0][0] = cosa * cosa * comm;
				stiffmatrix[0][1] = cosa * sina * comm;
				stiffmatrix[0][2] = -stiffmatrix[0][0];
				stiffmatrix[0][3] = -stiffmatrix[0][1];
				stiffmatrix[1][0] =  stiffmatrix[0][1];
				stiffmatrix[1][1] = sina * sina * comm;
				stiffmatrix[1][2] = -stiffmatrix[0][1];
				stiffmatrix[1][3] = -stiffmatrix[1][1];
				stiffmatrix[2][0] =  stiffmatrix[0][2];
				stiffmatrix[2][1] =  stiffmatrix[1][2];
				stiffmatrix[2][2] =  stiffmatrix[0][0];
				stiffmatrix[2][3] =  stiffmatrix[0][1];
				stiffmatrix[3][0] =  stiffmatrix[0][3];
				stiffmatrix[3][1] =  stiffmatrix[1][3];
				stiffmatrix[3][2] =  stiffmatrix[2][3];
				stiffmatrix[3][3] =  stiffmatrix[1][1];

				/* Assembly of Element Stiffness to Overall
				   Stiffness   */
				lm[1] = 2 * member[m][0] - 1;
				lm[0] = lm[1] - 1;
				lm[3] = 2 * member[m][1] - 1;
				lm[2] = lm[3] - 1;
				for (k = 0; k < 4; k++){
					for(g=0; g<4; g++){
						i = lm [k];
						j = lm [g];
						a[i][j] = a[i][j] + stiffmatrix[k][g];
					}
				}

			}
			/*end stiff*/

			/*displ*/
			
			for(n = 0; n < num_node; n++){
				nn = (2 * n) + 1;
				b[nn] = fY[n];
				b[nn-1] = fX[n];
			}

			for(n = 0; n < num_node; n++){
				nn = (2 * n) + 1;
				switch(kode[n]){
					case 1:
						for(i=0;i<neq;i++){
							b[i]=b[i]-(a[i][nn-1]*fX[n]);
							a[i][nn-1]=0.;
							a[nn-1][i]=0.;
						}
						a[nn-1][nn-1]=1.00;
						b[nn-1]=fX[n];
						break;
					case 2:
						for(i=0;i<neq;i++){
							b[i]=b[i]-(a[i][nn]*fY[n]);
							a[i][nn]=0.;
							a[nn][i]=0.;
						}
						a[nn][nn]=1.00;
						b[nn]=fY[n];
						break;
					case 3:
						for(i=0;i<neq;i++){
							b[i]=b[i]-a[i][nn-1]*fX[n];
							a[i][nn-1]=0.;
							a[nn-1][i]=0.;
						}
						a[nn-1][nn-1]=1.00;
						b[nn-1]=fX[n];
						for(i=0;i<neq;i++){
							b[i]=b[i]-a[i][nn]*fY[n];
							a[i][nn]=0.;
							a[nn][i]=0.;
						}
						a[nn][nn]=1.00;
						b[nn]=fY[n];
						break;
				}
			}
			
			/*	symsol();*/
			neq1 = nl = 0;
			neq1=neq-1;
			nl = neq - 1;

			for (n=0;n<nl;n++){
				if (a[n][n]<=0.){
					printf("zero or negative main-diagonal %d\n",n);
					exit(2);
				}

				n1 = n +1 ;
				for (j=n1;j<neq;j++) {
					a[n][j]=a[n][j]/a[n][n];
				}
				for(i=n1;i<neq;i++)
				{
					if(a[n][i]==0.0) b[i]=b[i]-a[n][i]*b[n];
					else{
						for(j=i;j<neq;j++){
							a[i][j]=a[i][j]-a[i][n]*a[n][j];
							a[j][i]=a[i][j];
						}
					}
					b[i]=b[i]-a[n][i]*b[n];
				}
				b[n]=b[n]/a[n][n];
			}
			
			m=neq1;
			b[m]=b[m]/a[m][m];
			for(n=0;n<nl;n++){
				m1=m;
				m=m-1;
				for(j=m1;j<neq;j++) {
					b[m]=b[m]-b[j]*a[m][j];
				}
			}
			/*end symbol*/
			
			for(i=0;i<num_node;i++){
				displX[i] = b[2*i];
				displY[i] = b[2*i+1];
			}
			/*end displ*/


			/*stress*/
			for(m=0;m<num_elem;m++){
				i   =  member[m][0] - 1;
				j   =  member[m][1] - 1;
				dx  =  coordx[j]-coordx[i];
				dy  =  coordy[j]-coordy[i];
				xl  =  sqrt((dx*dx) + (dy * dy));
				du  =  displX[j] -displX[i];
				dv  =  displY[j] -displY[i];
				dl  =  dv*dy/xl+du*dx/xl;
				stress[m] =  dl*ymod[m]/xl;
				p   =  area[m]*stress[m];
			}
			/*end stress*/
			
			/*stress constraint (25 ksi)*/
			for(i = 0; i < 10; i++) {
				y[i] = fabs(stress[i]);
				if(y[i] > 25){
					y[i] -=25;
					
				}else{
					y[i] = 0;
				}
				pen += y[i];      
			}
			/*displacement constraint (2 in.)*/
			for(i = 0; i < 6; i++) {
				y[10+i] = fabs(displX[i])/* * 1000.*/;
				if(y[10+i] > 2) y[10+i] -= 2;
				else y[10+i] = 0;
				pen += y[10+i];
			}

			/*displacement constraint (2 in.)*/
			for(i = 0; i < 6; i++) {
				y[16+i] = fabs(displY[i])/* * 1000.*/;
				if(y[16+i] > 2) y[16+i] -= 2;
				else y[16+i] = 0;
				pen += y[16+i];
			}
			
			pen *=P;
			free(area);
			break;
	}
	if(cond) memcpy(c_f,y,rest*sizeof(double));
	return pen;
}

double objfunc(double *sol, int cond){

    int j, i;
    double top = 0.00 , top1 = 0.00, top2 = 0.00;
    double aux = 0.0;
    double aux1 = 0.0;
    double somF = 0.;
    double ro = 0.1;
	double *area;
	double length[10] = {
	360.,
	360.,
	360.,
	360.,
	360.,
	360.,
	509.116882454,
	509.116882454,
	509.116882454,
	509.116882454
	};
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

    switch (FUNCTION) {
    case 0: //Rastrigin

        for(j=0;j<DIM;j++){
            top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
        }
        return top;

    case 1: //Schaffer

        top1 = 0;
        for(j=0;j<DIM;j++){
        	top=top+(pow(sol[j],(double)2));
        }
        top = pow(top,(double)0.25);
        for(j=0;j<DIM;j++){
        	top1=top1+(pow(sol[j],(double)2));
        }
        top1=pow(top1,(double)0.1);
        top1 = pow(sin(50*top1),(double)2) +1.0;

        return top*top1;

    case 2: //Griewank

        top=0;
        top1=0;
        top2=1;
        for(j=0;j<DIM;j++){
        	top1=top1+pow((sol[j]),(double)2);
        	top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);
        }
        top=(1/(double)4000)*top1-top2+1;

        return top;

    case 3: //Ackley

        for (i = 0; i < DIM; i++){
        	aux += sol[i]*sol[i];
        }
        for (i = 0; i < DIM; i++){
        	aux1 += cos(2.0*M_PI*sol[i]);
        }

        return (-20.0*(exp(-0.2*sqrt(1.0/(float)DIM*aux)))-exp(1.0/(float)DIM*aux1)+20.0+exp(1));

    case 4: //Rosenbrock

        for (i = 0; i < DIM-1; i++){
            top=top+100.*pow((sol[i+1] - pow(sol[i],2.)),2) + pow((sol[i] - 1.),2);
        }
        

       return top;

    case 5: //Sphere
		for(j=0;j<DIM;j++){
			top=top+sol[j]*sol[j];
		}
		return top;

    case 6: //Michaelewicz10
		for(j=0;j<DIM;j++){
			top=top+sin(sol[j])*pow(sin(((j+1)*(sol[j]*sol[j]))/M_PI),20);
		}
		return -top;

    case 7: //Booth
		return pow((sol[0]+(2*sol[1])-7),2)+pow((2*sol[0]+sol[1]-5),2);

    case 8: //Quartic
 		for(j=0;j<DIM;j++){
			top+=((j+1)*pow(sol[j],4));
		}
		return top+randon(0,1);
	
	case 9:
		top=constr(sol,cond);
		return (0.0624*(sol[0]+sol[1]+sol[2]+sol[3]+sol[4]))+top;
	
	case 10:
		top=constr(sol,cond);
		return (5000/(((sol[2]*(sol[1]-2*sol[3]))/12)+((sol[0]*pow(sol[3],3))/6)+(2*sol[0]*sol[3]*0.25*pow(sol[1]-sol[3],2))))+top;
	
	case 11:
		top=constr(sol,cond);
		if(cond) printf("Penalty: %g\n", top);
		return (1.10471*pow(sol[0],2)*sol[1]+0.04811*sol[2]*sol[3]*(14.0+sol[1]))+top;
	
	case 12://Pressure Vessel
		top=constr(sol,cond);
		if(cond) printf("Penalty: %g\n", top);
		return ((0.6224*sol[0]*sol[2]*sol[3])+(1.7781*sol[1]*pow(sol[2],2))+(3.1661*pow(sol[0],2)*sol[3])+(19.84*pow(sol[0],2)*sol[2]))+top;
	
	case 13: //Tension/compression string
		top=constr(sol,cond);
		return ((sol[2]+2)*sol[1]*pow(sol[0],2))+top;
	
	case 14: //Speed Reducer
		top=constr(sol,cond);
		return (0.7854*sol[0]*pow(sol[1],2)*(3.3333*pow(sol[2],2)+14.9334*sol[2]-43.0934)-1.508*sol[0]*(pow(sol[5],2)+pow(sol[6],2))+7.477*(pow(sol[5],3)+pow(sol[6],3))+0.7854*(sol[3]*pow(sol[5],2)+sol[4]*pow(sol[6],2)))+top;
    case 15: 		
		top=constr(sol,cond);		    	
    		area = (double*) malloc (DIM * sizeof(double));
		if(area==NULL)exit(0);
		for (i=0; i < DIM; i++){
			area[i] = sol[i];
		}
		/*weight of the truss*/
		for (i = 0; i < 10; i++) {
			somF += ro*area[i]*length[i];
		}
		
		return somF+top;
	default:
        printf("Info: Invalid function..\n") ;
        exit(0);
    }
}

void *th_init_pop(void *argThread){
	int j,k;
	slice *s = (slice*)argThread;
	
	if(FUNCTION!=10 && FUNCTION!=11 && FUNCTION!=12 && FUNCTION!=13 && FUNCTION!=14){
		for (j=s->inicio;j<s->fim;j++){//each individual
			for (k=0; k<DIM;k++){ //each dimension of the individual
				pop[j][k] = randon(lb[0],ub[0]);
			}	
			fo[j] = objfunc(pop[j], 0);
		}
	}else{
		for (j=s->inicio;j<s->fim;j++){//each individual
			for (k=0; k<DIM;k++){ //each dimension of the individual
				pop[j][k] = randon(lb[k],ub[k]);
			}
			fo[j] = objfunc(pop[j], 0);
		}
	}
}

void *th_sos(void* argThread){
	int i,k;
	slice *s = (slice*)argThread;

	for(i=s->inicio;i<s->fim;i++){
		mutualism_phase(i,pop,best,fo);
		num_fit_eval+=2;
		commensalism_phase(i,pop,best,fo);
		num_fit_eval++;
		parasitism_phase(i,pop,fo);
		num_fit_eval++;	
			
		for(k=0;k<POP_SIZE;k++){
			if(fo[k] <= bestfo){
			        pthread_mutex_lock(&data_mutex);
                                bestfo=fo[k];
				best_index=k;
                                pthread_mutex_unlock(&data_mutex);
			}
		}
		for(k=0;k<DIM;k++){
                    pthread_mutex_lock(&data_mutex);
                    best[k]=pop[best_index][k];
                    pthread_mutex_unlock(&data_mutex);
                }
	}
}

void showConst(double *var, int r, FILE *file){
	int i;
	if(COND){
		if(constr(best,1)==0)var[r]=bestfo;
		else var[r]=2147483646;
		//values of constraints
		
		switch(FUNCTION){
			case 9: //Cantilever Beam
				fprintf(file,"g1=%g ",c_f[0]);
				if(c_f[0]>1) fprintf(file, "Fail\n");
				else fprintf(file, "Ok\n");
				printf("g1=%g ",c_f[0]);
				if(c_f[0]>1)printf("Fail\n");
				else printf("Ok\n");
				break;
			case 10: //I-Beam vertical deflection 
				fprintf(file, "g1=%g ",c_f[0]);
				if(c_f[0]>300) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file, "g2=%g ",c_f[1]);
				if(c_f[1]>56) fprintf(file, "Fail\n");
				else fprintf(file, "Ok\n");
				printf("g1=%g ",c_f[0]);
				if(c_f[0]>300) printf("Fail ");
				else printf("Ok ");
				printf("g2=%g ",c_f[1]);
				if(c_f[1]>56) printf("Fail\n");
				else printf("Ok\n");
				break;
			case 11: //Welded Beam 
				fprintf(file,"g1=%g ",c_f[0]);fprintf(file,"g1=%g ",c_f[0]);
				if(c_f[0]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g2=%g ",c_f[1]);
				if(c_f[1]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g3=%g ",c_f[2]);
				if(c_f[2]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g4=%g ",c_f[3]);
				if(c_f[3]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g5=%g ",c_f[4]);
				if(c_f[4]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g6=%g ",c_f[5]);
				if(c_f[5]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g7=%g ",c_f[6]);
				if(c_f[6]>0) fprintf(file, "Fail\n");
				else fprintf(file, "Ok\n");
				printf("g1=%g ",c_f[0]);
				if(c_f[0]>0) printf("Fail ");
				else printf("Ok ");
				printf("g2=%g ",c_f[1]);
				if(c_f[1]>0) printf("Fail ");
				else printf("Ok ");
				printf("g3=%g ",c_f[2]);
				if(c_f[2]>0) printf("Fail ");
				else printf("Ok ");
				printf("g4=%g ",c_f[3]);
				if(c_f[3]>0) printf("Fail ");
				else printf("Ok ");
				printf("g5=%g ",c_f[4]);
				if(c_f[4]>0) printf("Fail ");
				else printf("Ok ");
				printf("g6=%g ",c_f[5]);
				if(c_f[5]>0) printf("Fail ");
				else printf("Ok ");
				printf("g7=%g ",c_f[6]);
				if(c_f[6]>0) printf("Fail\n");
				else printf("Ok\n");
				break;
			case 12: //Pressure Vessel 
				fprintf(file,"g1=%g ",c_f[0]);
				if(c_f[0]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g2=%g ",c_f[1]);
				if(c_f[1]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g3=%g ",c_f[2]);
				if(c_f[2]>0) fprintf(file, "Fail ");
				else fprintf(file, "Ok ");
				fprintf(file,"g4=%g ",c_f[3]);
				if(c_f[3]>0) fprintf(file, "Fail\n");
				else fprintf(file, "Ok\n");
				printf("g1=%g ",c_f[0]);
				if(c_f[0]>0) printf("Fail ");
				else printf("Ok ");
				printf("g2=%g ",c_f[1]);
				if(c_f[1]>0) printf("Fail ");
				else printf("Ok ");
				printf("g3=%g ",c_f[2]);
				if(c_f[2]>0) printf("Fail ");
				else printf("Ok ");
				printf("g4=%g ",c_f[3]);
				if(c_f[3]>0) printf("Fail\n");
				else printf("Ok\n");
				break;
			case 13: //Tension/compression string
				fprintf(file,"g1=%g g2=%g g3=%g g4=%g\n",c_f[0],c_f[1],c_f[3],c_f[3]);
				printf("g1=%g g2=%g g3=%g g4=%g\n",c_f[0],c_f[1],c_f[3],c_f[3]);
				break;
			case 14: //Speed Reducer(Gear Train)
				fprintf(file,"g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",c_f[0],c_f[1],c_f[3],c_f[3],c_f[4],c_f[5],c_f[6],c_f[7],c_f[8],c_f[9],c_f[10]);
				printf("g1=%g g2=%g g3=%g g4=%g g5=%g g6=%g g7=%g g8=%g g9=%g g10=%g g11=%g\n",c_f[0],c_f[1],c_f[3],c_f[3],c_f[4],c_f[5],c_f[6],c_f[7],c_f[8],c_f[9],c_f[10]);
				break;
			case 15: //10-Bar-Truss
				printf("Stress Violation: ");
				for(i=0;i<10;i++){
					printf("g(%i)=%g ",i+1, c_f[i]);
				}
				printf("\n");
				printf("DispX: ");
				for(i=10;i<16;i++){
					printf("g(%i)=%g ",i+1, c_f[i]);
				}
				printf("\n");
				printf("DispY: ");
				for(i=16;i<22;i++){
					printf("g(%i)=%g ",i+1,c_f[i]);
				}
				printf("\n");

				fprintf(file,"Stress Violation: ");
				for(i=0;i<10;i++){
					fprintf(file,"g(%i)=%g ",i+1, c_f[i]);
				}
				fprintf(file,"\n");
				fprintf(file,"DispX: ");
				for(i=10;i<16;i++){
					fprintf(file,"g(%i)=%g ",i+1, c_f[i]);
				}
				fprintf(file,"\n");
				fprintf(file,"DispY: ");
				for(i=16;i<22;i++){
					fprintf(file,"g(%i)=%g ",i+1,c_f[i]);
				}
				fprintf(file,"\n");
				break;
		}
	}
}

void sos_iter(){
	int i;	
	pthread_t *threads;

	threads = (pthread_t*)malloc(CORES*sizeof(pthread_t));
	if(threads==NULL)exit(0);
	
	for (i = 0; i < CORES; ++i){
		pthread_create(&threads[i], NULL, th_sos, (void *) (argThread+i));
	}
		
	for (i = 0; i < CORES; i++){
      	pthread_join(threads[i], NULL);
   	}
   	free(threads);
}


void free_slice(){
	free(argThread);	
}


void initPop(){
	int i;	
	int pedaco = ((double)POP_SIZE/(double)CORES); 
	pthread_t *threads;
	
	threads = (pthread_t*)malloc(CORES*sizeof(pthread_t));
	if(threads==NULL)exit(0);
	
	for(i=0;i<CORES;i++){
		argThread[i].tid = i;
		argThread[i].inicio = i*pedaco;
		if(i==CORES-1){
			argThread[i].fim = POP_SIZE;
		}else argThread[i].fim = argThread[i].inicio+(pedaco);
		pthread_create(&threads[i], NULL, th_init_pop, (void *) (argThread+i));
	}
			
	for (i = 0; i < CORES; i++){
    		pthread_join(threads[i], NULL);
   	}

   		
   	free(threads);
}

void mutualism_phase(int index_i, double **pop_th, double *best_th, double *fo_th){
	int i;
	int bf1=0,bf2=0;	
	double *mutual;
	double *new_x_i,*new_x_j;
	double newfo_i,newfo_j;
	double *array_rand;
        double *l,*u;
	int index_j=index_i;

	while(index_j==index_i){
		index_j=(int)randon(0,POP_SIZE);
	}//pick xj	

        l = (double*)malloc(DIM*sizeof(double));
        u = (double*)malloc(DIM*sizeof(double));
        memcpy(l,lb,DIM*sizeof(double));
        memcpy(u,ub,DIM*sizeof(double));
	array_rand=(double*)malloc(DIM*sizeof(double));
	new_x_i = (double*)malloc(DIM*sizeof(double));
	new_x_j = (double*)malloc(DIM*sizeof(double));	
	mutual = (double*)malloc(DIM*sizeof(double));

	
	bf1=(int)randon(1,2.5);//benefit factor1
	bf2=(int)randon(1,2.5);//benefit factor2
	for(i=0;i<DIM;i++){
		array_rand[i]=randon(0,1);
	}
	
	//cáculo Eqs(1) (2) (3)
	for(i=0;i<DIM;i++){
		pthread_mutex_lock(&data_mutex);
                mutual[i]=(pop_th[index_i][i]+pop_th[index_j][i])/2;
		new_x_i[i]=pop_th[index_i][i]+(array_rand[i]*(best_th[i]-(mutual[i]*bf1)));
		new_x_j[i]=pop_th[index_j][i]+(array_rand[i]*(best_th[i]-(mutual[i]*bf2)));
		pthread_mutex_unlock(&data_mutex);
                if(FUNCTION!=10 && FUNCTION!=11 && FUNCTION!=12 && FUNCTION!=13 && FUNCTION!=14){
			if(new_x_i[i]>u[0])new_x_i[i]=u[0];
			if(new_x_j[i]>u[0])new_x_j[i]=u[0];
			if(new_x_i[i]<l[0])new_x_i[i]=l[0];
			if(new_x_j[i]<l[0])new_x_j[i]=l[0];
		}else{
			if(new_x_i[i]>u[i])new_x_i[i]=u[i];
			if(new_x_j[i]>u[i])new_x_j[i]=u[i];
			if(new_x_i[i]<l[i])new_x_i[i]=l[i];
			if(new_x_j[i]<l[i])new_x_j[i]=l[i];
		}
	}
	//	//
	newfo_i=objfunc(new_x_i, 0);//calculates fitness for new_x_i
	newfo_j=objfunc(new_x_j, 0);//calculates fitness for new_x_j
	pthread_mutex_lock(&data_mutex);
	if(fo_th[index_i]>=newfo_i){//greedy selection for xi
		for(i=0;i<DIM;i++){
			pop_th[index_i][i]=new_x_i[i];
		}
		fo_th[index_i] = newfo_i;
	}
	pthread_mutex_unlock(&data_mutex);
	pthread_mutex_lock(&data_mutex);
	if(fo_th[index_j]>=newfo_j){//greedy selection for xj
		for(i=0;i<DIM;i++){
			pop_th[index_j][i]=new_x_j[i];
		}
		fo_th[index_j]=newfo_j;
	}
	pthread_mutex_unlock(&data_mutex);
	free(array_rand);
	free(new_x_i);	
	free(new_x_j);
	free(mutual);		
}

void commensalism_phase(int index_i, double **pop_th, double *best_th, double *fo_th){
	int i;
	double *array_rand;
	int index_j=index_i;
	double *new_x_i,*u,*l;
	double newfo_i;

	while(index_j==index_i){
		index_j=(int)randon(0, POP_SIZE);
	}//pick xj	

	//alloc arrays//
        l = (double*)malloc(DIM*sizeof(double));
        u = (double*)malloc(DIM*sizeof(double));
	array_rand=(double*)malloc(DIM*sizeof(double));
	new_x_i = (double*)malloc(DIM*sizeof(double));
	//
        memcpy(l,lb,DIM*sizeof(double));
        memcpy(u,ub,DIM*sizeof(double));

	//puts values for array_rand and new_x_i
	for(i=0;i<DIM;i++){
		array_rand[i]=randon(-1,1);
		//pthread_mutex_lock(&data_mutex);
                new_x_i[i]=pop_th[index_i][i]+(array_rand[i]*(best_th[i]-pop_th[index_j][i]));
		//pthread_mutex_unlock(&data_mutex);
                if(FUNCTION!=10 && FUNCTION!=11 && FUNCTION!=12 && FUNCTION!=13 && FUNCTION!=14){
			if(new_x_i[i]<l[0])new_x_i[i]=l[0];
			if(new_x_i[i]>u[0])new_x_i[i]=u[0];
		}else{
			if(new_x_i[i]<l[i])new_x_i[i]=l[i];
			if(new_x_i[i]>u[i])new_x_i[i]=u[i];
		}
	}
	

	newfo_i=objfunc(new_x_i, 0);//calculates fitness for new_x_i
	
	
	if(fo_th[index_i]>=newfo_i){////greedy selection for xi
		for(i=0;i<DIM;i++){
			//pthread_mutex_lock(&data_mutex);
                        pop_th[index_i][i]=new_x_i[i];
		        //pthread_mutex_unlock(&data_mutex);
                }
			//pthread_mutex_lock(&data_mutex);
                        fo_th[index_i]=newfo_i;
	                //pthread_mutex_unlock(&data_mutex);                
        }

	free(array_rand);
	free(new_x_i);	
}

void parasitism_phase(int index_i, double **pop_th, double *fo_th){
	double *parasite,*u,*l;
	double parasite_fo;
	int pick;
	int index_j=index_i;
	int i;

	parasite=(double*)malloc(DIM*sizeof(double));//alloc parasite array	
	pick=(int)randon(0,DIM);//chooses the DIMension for change
	l = (double*)malloc(DIM*sizeof(double));
        u = (double*)malloc(DIM*sizeof(double));
        memcpy(l,lb,DIM*sizeof(double));
        memcpy(u,ub,DIM*sizeof(double));

	
	for(i=0;i<DIM;i++){//copies xi for parasite
            //pthread_mutex_lock(&data_mutex);		
            parasite[i]=pop_th[index_i][i];
	    //pthread_mutex_unlock(&data_mutex);
        }
	

	while(index_j==index_i){
		index_j=(int)randon(0,POP_SIZE);
	}
	
	if(FUNCTION!=10 && FUNCTION!=11 && FUNCTION!=12 && FUNCTION!=13 && FUNCTION!=14){
		parasite[pick]=randon(l[0],u[0]);//change the value of the DIMension choosen
	}else{
		parasite[pick]=randon(l[pick],u[pick]);//change the value of the DIMension choosen
	}
	parasite_fo=objfunc(parasite, 0);//calculates fitness for parasite
	
	if(fo_th[index_j]>=parasite_fo){//greedy selection between xj and parasite
		for(i=0;i<DIM;i++){
		    //pthread_mutex_lock(&data_mutex);	
                    pop_th[index_j][i]=parasite[i];
		    //pthread_mutex_unlock(&data_mutex);
                }
		//pthread_mutex_lock(&data_mutex); 
                fo_th[index_j]=parasite_fo;
		//pthread_mutex_unlock(&data_mutex);		
        }
	free(parasite);	
}

int AvgStdDev(double *Avg,double *StdDev,double Var[]){
	int i;
	int qtd=0;

	*Avg = 0;
	*StdDev = 0;
	for (i=0;i<RUN;i++){
		if(Var[i]!=2147483646){
			*Avg += Var[i];
			qtd++;
		}
	}
	*Avg /= qtd;

	for (i=0;i<RUN;i++){
		if(Var[i]!=2147483646){
			*StdDev += pow((*Avg-Var[i]),2);
		}
	}
	*StdDev /= qtd;
   *StdDev = sqrt(*StdDev);
   return qtd;
}

void plot(FILE *shellComands, char *run){
	if(strcmp(run,"Final")==0){
		shellComands = popen ("gnuplot -persistent", "w");
		fprintf(shellComands,"	reset\n"); 
		fprintf(shellComands,"	set terminal pngcairo\n");
		fprintf(shellComands,"	set output \"dadosplot/sosplotMedia_M-Media_Bfo-%s.png\"\n",run);
		fprintf(shellComands,"	set grid\n");
		fprintf(shellComands,"	set title \"SOS Analysis\"\n");
		fprintf(shellComands,"	set xlabel \"Itera\303\247\303\265es\"\n");
		fprintf(shellComands,"	set ylabel \"Fo\" \n");
		fprintf(shellComands,"	plot \"dadosplot//dadosplot%s.txt\" using 1:3 title \"Media_M\" w lines, \"dadosplot//dadosplot%s.txt\" using 1:2 title \"Media_Bfo\" w lines",run,run);
		pclose(shellComands); 
	}else{
		shellComands = popen ("gnuplot -persistent", "w");
		fprintf(shellComands,"	reset\n"); 
		fprintf(shellComands,"	set terminal pngcairo\n");
		fprintf(shellComands,"	set output \"dadosplot/sosplotMedia_fo-Best_Fo-%s.png\"\n",run);
		fprintf(shellComands,"	set grid\n");
		fprintf(shellComands,"	set title \"SOS Analysis\"\n");
		fprintf(shellComands,"	set xlabel \"Itera\303\247\303\265es\"\n");
		fprintf(shellComands,"	set ylabel \"Fo\" \n");
		fprintf(shellComands,"	plot \"dadosplot//dadosplot%s.txt\" using 1:3 title \"Media Fo\" w lines , \"dadosplot//dadosplot%s.txt\" using 1:2 title \"best_thFo\" w lines",run,run);
		pclose(shellComands); 	
	}
}

char *converteDecChar(char *strf, int dec){
	int i;
	char str[100];
	if(dec==0){
		strf[0]='0';
		strf[1]='\0';
		return strf;
	}
	for ( i = 0; i < 100; i++ ){
		if ( dec == 0 ) break;
		str[i] = (dec % 10) + 48;

		dec /= 10;
	}
	str[i]='\0';
	for ( i = strlen(str)-1; i >= 0; i-- ){
		strf[(strlen(str)-1)-i]=str[i];
	}
	strf[strlen(str)]='\0';
	return strf;
}

int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer){
	char buffer[100];
	int len;
    int commentflag = 0;
	char *pch;
	char *pch2;
	
	do
	{
		if(fgets(buffer, 99, fp) == 0) return FAIL;
		buffer[99] = '\0';
		len = strlen(buffer);
	    if (buffer[len - 1] == '\n') buffer[len - 1] = '\0';

	    switch (commentflag) 
		{
		    case 0:
				if (strstr(buffer, "/*") != 0) 
				{
			    	commentflag = 1;
			    	if (strstr(buffer, "*/") != 0)
						commentflag = 2;
				}
				break;

		    case 1:
				if (strstr(buffer, "*/") != 0)
				    commentflag = 2;
				break;

			    case 2:
				if (strstr(buffer, "/*") != 0) 
				{
				    commentflag = 1;
				    if (strstr(buffer, "*/") != 0)
					commentflag = 2;
				}
				else
				    commentflag = 0;
				break;
			    }	
	}while(commentflag != 0);	

	//separate field name: token = "="
	if (strstr (buffer, fieldname) != 0)
	{
		pch = strtok (buffer,"=");

		while (pch != NULL)
		{
			pch2 = pch;
    		pch = strtok (NULL, "= ");
		}
		sscanf(pch2, format, inbuffer);
		return 1;//ok
	}
	else return 0; 
}

void AllocArrays(){
	int j;
	//int constn=0;

	argThread = (slice*)malloc(CORES*sizeof(slice));
	if(argThread==NULL)exit(0);

	c_f = (double*)malloc(rest*sizeof(double));
	if(c_f==NULL)exit(0);

	pop = (double**)malloc(POP_SIZE*sizeof(double*));
		double *dados = (double*)malloc(sizeof(double)*DIM*POP_SIZE);	
        for (j = 0;j < POP_SIZE; j++)
			pop[j] = &dados[j*DIM];


	fo = (double*)malloc (POP_SIZE*sizeof(double));

	best = (double*)malloc (DIM*sizeof(double));

	ub=(double*)malloc (DIM*sizeof(double));
	lb=(double*)malloc (DIM*sizeof(double));
}
