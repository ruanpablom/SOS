#include <stdio.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "sosfunctions.h"
#include "mersenne.h"

#define FAIL 0
//#define COND func!=10 && func!=11 && func!=12 && func!=13 && func!=14

double randon( double inferior, double superior){
	//  double aux = (float)inferior + ((superior - inferior)*rand()/(RAND_MAX+1.0));
	double aux = (double)inferior + ((superior - inferior)*MT_randInt(RAND_MAX)/(RAND_MAX+1.0));
	//double aux = (double)inferior + ((superior - inferior)* MT_randfloat()/(RAND_MAX+1.0));
	return aux;
}


/*void AllocArrays(int pop_size, int dim, double **pop, double *fo, double *best, double *ub, double *lb){
	int j;

	pop = (double**)malloc (pop_size*sizeof(double*));
        for (j = 0;j < pop_size; j++)
		pop[j] = (double*)malloc (dim * sizeof(double));

	fo = (double*)malloc (pop_size*sizeof(double));

	best = (double*)malloc (dim*sizeof(double));

	ub=(double*)malloc (dim*sizeof(double));
	lb=(double*)malloc (dim*sizeof(double));
}*/

void freeArrays(int pop_size, double **pop, double *fo, double *best, double *ub, double *lb){
	int i;
	free(fo);
	free(best);
	free(lb);
	free(ub);
	for (i = 0; i < pop_size; i++){
		free(pop[i]);
	}
	free(pop);
}

int GetParameters(char **argv, int *run, int *max_iter, int *pop_size, int *dim, int *function){
	FILE *file = fopen( argv[1], "r" );
    
	if (file == 0){
    	printf( "Could not open ini file! Usage ./<exec> <file.in>\n" );
		return -1;
	}else{
		ffscanf("RUN", file, "%d", run);
		ffscanf("MAX_ITER", file, "%d", max_iter);
		ffscanf("POP_SIZE", file, "%d", pop_size);
		ffscanf("DIM", file, "%d", dim); 
		ffscanf("FUNCTION",file, "%d", function);
	return 1;
    }
    fclose(file);
}

void showParameters(int func, int run, int max_iter, int pop_size, int dim){
	printf("***PARAMETERS***\n");
	printf("RUNS = %d\n", run);
	printf("MAX_ITER = %d\n", max_iter);
	printf("POP_SIZE = %d\n", pop_size);
	printf("DIM = %d\n", dim);
	switch (func){
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

void prepararObjFunc(int func, double *l, double *u){
    switch (func){
    	case 0: //Rastrigin
        	l[0] = -5.12;
        	u[0] = 5.12;
        	break;
    	case 1: //Schaffer
        	l[0] = -100.00;
        	u[0] = 100.00;
        	break;
    	case 2: //Griewank
       		l[0] = -600.00;
        	u[0] = 600.00;
        	break;
    	case 3: //Ackley
        	l[0] = -32.00;
        	u[0] = 32.00;
        	break;
    	case 4: //Rosenbrock
        	l[0] = -30.00;
        	u[0] = 30.00;
        	break;
    	case 5: //Sphere
			l[0] = -100.00;
			u[0] = 100.00;
			break;
    	case 6: //Michaelewicz10
			l[0] = 0.00;
			u[0] = M_PI;
			break;
    	case 7: //Booth
			l[0] = -10.00;
			u[0] = 10.00;
			break;
    	case 8: //Quartic
			l[0] = -1.28;
			u[0] = 1.28;
			break;
		case 9: //Cantilever beam
			l[0] = 0.01;
			u[0] = 100.0;
			break;
		case 10: //10: I-Beam vertical deflection 
			//b domain
			l[0] = 10.0;
			u[0] = 50.0;
			//
			//h domain
			l[1] = 10.0;
			u[1] = 80.0;
			//
			//tw domain
			l[2] = 0.9;
			u[2] = 5.0;
			//
			//tf domain
			l[3] = 0.9;
			u[3] = 5.0;
			//
			break;
		case 11: // Welded Beam 
			l[0] = 0.1;
			u[0] = 2.0;
			l[1] = 0.1;
			u[1] = 10.0;
			l[2] = 0.1;
			u[2] = 10.0;
			l[3] = 0.1;
			u[3] = 2.0;
			//
			break;
		case 12: // Pressure Vessel
			l[0] = 1.0;
			u[0] = 99.0;
			l[1] = 1.0;
			u[1] = 99.0;
			l[2] = 10.0;
			u[2] = 200.0;
			l[3] = 10.0;
			u[3] = 200.0;
			//
			break;
		case 13: // Tension/compression string
			l[0] = 0.05;
			u[0] = 2.0;
			l[1] = 0.25;
			u[1] = 1.3;
			l[2] = 2.0;
			u[2] = 15.0;
			//
			break;
		case 14: // Speed Reducer(Gear Train)
			l[0] = 2.6;
			u[0] = 3.6;
			l[1] = 0.7;
			u[1] = 0.8;
			l[2] = 17;
			u[2] = 28;
			l[3] = 7.3;
			u[3] = 8.3;
			l[4] = 7.3;
			u[4] = 8.3;
			l[5] = 2.9;
			u[5] = 3.9;
			l[6] = 5.0;
			u[6] = 5.5;
			//
			break;
		case 15:
			l[0] = 0.10;
			u[0] = 35.00;
			break;
    	default:
    	    printf("Info: Invalid function\n") ;
    	    exit(0);
    	}
}	

double constr(double sol[], int func, double *y){//calculate penalization
	double pen=0;
	double r,m,j,s,d,pc,t1,t2,t,p,l,e,g,t_max,s_max,d_max;
	switch(func){
		case 9:
			y[0]=((61/pow(sol[0],3))+(37/pow(sol[1],3))+(19/pow(sol[2],3))+(7/pow(sol[3],3))+(1/pow(sol[4],3)));
			if(y[0]>1)pen+=y[0];
			break;
		case 10:
			//sol[0]=b
    		//sol[1]=h
    		//sol[2]=tw
    		//sol[3]=tf
			y[0]=2*sol[0]*sol[2]+sol[2]*(sol[1]-(2*sol[3]))+1;
			if(y[0]>300)pen+=y[0];
			y[1]=(((18*sol[1]*10000)/(sol[2]*pow(sol[1]-(2*sol[3]),3)+2*sol[0]*sol[2]*(4*pow(sol[3],2)+3*sol[1]*(sol[1]-2*sol[3]))))+((15*sol[0]*1000)/(((sol[1]-2*sol[3])*pow(sol[2],3))+2*sol[2]*pow(sol[0],3))));
			if(y[1]>56)pen+=y[1];
			break;
		case 11:
			p = 6000; 
			l = 14; 
			e = 30e+6; 
			g = 12e+6;
    		t_max = 13600; 
    		s_max = 30000; 
    		d_max = 0.25;
			j=2*sqrt(2)*sol[0]*sol[1]*((pow(sol[1],2)/12)+0.25*pow((sol[0]+sol[2]),2));
			r=sqrt(0.25*(pow(sol[1],2)+pow((sol[0]+sol[2]),2)));
			m=p*(14+(sol[1]/2));
			t2=(m*r)/j;
			t1=p/(sqrt(2)*sol[0]*sol[1]);
			t=sqrt(pow(t1,2)+2*t1*t2*(sol[1]/(2*r))+pow(t2,2));
			s=(6*p*l)/(sol[3]*pow(sol[2],2));
			d=(4*p*pow(l,3))/(sol[3]*pow(sol[2],3)*e);
			pc=((4.013*e*sqrt((pow(sol[2],2)*pow(sol[3],6))/36))/196)*(1-((sol[2]*sqrt(e/(4*g)))/28));
			y[0]=t-t_max;
			if(ceil(y[0])>0)pen+=y[0];
			y[1]=s-s_max;
			if(ceil(y[1])>0)pen+=y[1];
			y[2]=sol[0]-sol[3];
			if(ceil(y[2])>0)pen+=y[2];
			y[3]=1.10471*pow(sol[0],2)+0.04811*sol[2]*sol[3]*(14+sol[1])-5;
			if(ceil(y[3])>0)pen+=y[3];
			y[4]=0.125-sol[0];
			if(ceil(y[4])>0)pen+=y[4];
			y[5]=d-d_max;
			if(ceil(y[5])>0)pen+=y[5];
			y[6]=6000-pc;
			if(ceil(y[6])>0)pen+=y[6];
			break;
		case 12: //Pressure Vessel
			y[0]=-sol[0]+0.0193*sol[2];
			if(!(ceil(y[0])<=0))pen+=y[0];
			y[1]=-sol[1]+0.00954*sol[2];
			if(!(ceil(y[1])<=0))pen+=y[1];
			y[2]=-(M_PI)*pow(sol[2],2)*sol[3]-(4*pow(sol[2],3))*(M_PI)/3+1296000;
			if(!(ceil(y[2])<=0))pen+=y[2];
			y[3]=sol[3]-240;
			if(!(ceil(y[3])<=0))pen+=y[3];
		break;
		case 13: //Tension/compression string
			y[0]=1-((pow(sol[1],3)*sol[2])/(71785*pow(sol[0],4)));
			if(y[0]>0)pen+=y[0];
			y[1]=((((4*pow(sol[1],2))-(sol[0]*sol[1]))/(12566*((sol[1]*pow(sol[0],3))-pow(sol[0],4))))+(1/(5108*pow(sol[0],2)))-1);
			if(y[1]>0)pen+=y[1];
			y[2]=1-((140.45*sol[0])/(sol[2]*pow(sol[1],2)));
			if(y[2]>0)pen+=y[2];
			y[3]=((sol[0]+sol[1])/1.5)-1;
			if(y[3]>0)pen+=y[3];
			break;
		case 14: // Speed Reducer(Gear Train)
			y[0]=(27/(sol[0]*pow(sol[1],2)*sol[2]))-1;
			if(y[0]>0)pen+=y[0];//
			y[1]=(397.5/(sol[0]*pow(sol[1],2)*pow(sol[2],2)))-1;
			if(y[1]>0)pen+=y[1];
			y[2]=((1.93*pow(sol[3],3))/(sol[1]*sol[2]*pow(sol[5],4)))-1;
			if(y[2]>0)pen+=y[2];
			y[3]=((1.93*pow(sol[4],3))/(sol[1]*sol[2]*pow(sol[6],4)))-1;
			if(y[3]>0)pen+=y[3];
			y[4]=(sqrt(pow(((745*sol[3])/(sol[1]*sol[2])),2)+(16.9*pow(10,6))))*(1/(110*pow(sol[5],3)))-1;
			if(y[4]>0)pen+=y[4];
			y[5]=(sqrt(pow(((745*sol[4])/(sol[1]*sol[2])),2)+(157.5*pow(10,6))))*(1/(85*pow(sol[6],3)))-1;
			if(y[5]>0)pen+=y[5];
			y[6]=((sol[1]*sol[2])/40)-1;
			if(y[6]>0)pen+=y[6];
			y[7]=((5*sol[1])/sol[0])-1;
			if(y[7]>0)pen+=y[7];
			y[8]=(sol[0]/(sol[1]*12))-1;
			if(y[8]>0)pen+=y[8];
			y[9]=((1.5*sol[5]+1.9)/sol[3])-1;
			if(y[9]>0)pen+=y[9];
			y[10]=((1.1*sol[6]+1.9)/sol[4])-1;
			if(y[10]>0)pen+=y[10];
			break;
	}
	return pen;
}

double objfunc(double sol[], double y[], int func, int dim){
    int j, i;
    double top = 0.00 , top1 = 0.00, top2 = 0.00;
    double aux = 0.0;
    double aux1 = 0.0;

    //for 10-bar-truss
	double *area;
	double stress[10] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
	double displX[6] = {0., 0., 0., 0., 0., 0.};
	double displY[6] = {0., 0., 0., 0., 0., 0.};
	double somF = 0.;
	double somg1 = 0.;
	double somg2 = 0.;

	double ro = 0.1;
	double P  = 10000.;

	double E  = 10000.;

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
	double du = 0.,dv = 0.,p = 0.,dl = 0.;

	int nl = 0;
	int n1 = 0,m1 = 0;
	int neq1 = 0;

	int m = 0, k = 0, g = 0;
	double cosa = 0.,sina = 0.,comm = 0.;

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

    switch (func) {
    case 0: //Rastrigin

        for(j=0;j<dim;j++)
        {
            top=top+(pow(sol[j],(double)2)-10*cos(2*M_PI*sol[j])+10);
        }
        return top;

    case 1: //Schaffer

        top1 = 0;
        for(j=0;j<dim;j++)
        {
        top=top+(pow(sol[j],(double)2));
        }
        top = pow(top,(double)0.25);
        for(j=0;j<dim;j++)
        {
        top1=top1+(pow(sol[j],(double)2));
        }
        top1=pow(top1,(double)0.1);
        top1 = pow(sin(50*top1),(double)2) +1.0;

        return top*top1;

    case 2: //Griewank

        top=0;
        top1=0;
        top2=1;
        for(j=0;j<dim;j++)
        {
        top1=top1+pow((sol[j]),(double)2);
        top2=top2*cos((((sol[j])/sqrt((double)(j+1)))*M_PI)/180);
        }
        top=(1/(double)4000)*top1-top2+1;

        return top;

    case 3: //Ackley

        for (i = 0; i < dim; i++)
        {
        aux += sol[i]*sol[i];
        }
        for (i = 0; i < dim; i++)
        {
        aux1 += cos(2.0*M_PI*sol[i]);
        }

        return (-20.0*(exp(-0.2*sqrt(1.0/(float)dim*aux)))-exp(1.0/(float)dim*aux1)+20.0+exp(1));

    case 4: //Rosenbrock

        for (i = 0; i < dim-1; i++)
        {
            top=top+100.*pow((sol[i+1] - pow(sol[i],2.)),2) + pow((sol[i] - 1.),2);
        }
        

       return top;

    case 5: //Sphere
	for(j=0;j<dim;j++)
	{
		top=top+sol[j]*sol[j];
	}

	return top;
    case 6: //Michaelewicz10
	for(j=0;j<dim;j++)
	{
		top=top+sin(sol[j])*pow(sin(((j+1)*(sol[j]*sol[j]))/M_PI),20);
	}

	return -top;
    case 7: //Booth
	return pow((sol[0]+(2*sol[1])-7),2)+pow((2*sol[0]+sol[1]-5),2);

    case 8: //Quartic
 		for(j=0;j<dim;j++){
			top+=((j+1)*pow(sol[j],4));
		}
		return top+randon(0,1);
	case 9:
		top=constr(sol, func,y);
		return (0.0624*(sol[0]+sol[1]+sol[2]+sol[3]+sol[4]))+top;
	case 10:
		top=constr(sol, func,y);
		return (5000/(((sol[2]*(sol[1]-2*sol[3]))/12)+((sol[0]*pow(sol[3],3))/6)+(2*sol[0]*sol[3]*0.25*pow(sol[1]-sol[3],2))))+top;
	case 11:
		top=constr(sol, func,y);
		return (1.10471*pow(sol[0],2)*sol[1]+0.04811*sol[2]*sol[3]*(14.0+sol[1]))+top;
	case 12://Pressure Vessel
		top=constr(sol, func,y);
		return ((0.6224*sol[0]*sol[2]*sol[3])+(1.7781*sol[1]*pow(sol[2],2))+(3.1661*pow(sol[0],2)*sol[3])+(19.84*pow(sol[0],2)*sol[2]))+top;
	case 13: //Tension/compression string
		top=constr(sol, func,y);
		return ((sol[2]+2)*sol[1]*pow(sol[0],2))+top;
	case 14: //Tension/compression string
		top=constr(sol, func,y);
		return (0.7854*sol[0]*pow(sol[1],2)*(3.3333*pow(sol[2],2)+14.9334*sol[2]-43.0934)-1.508*sol[0]*(pow(sol[5],2)+pow(sol[6],2))+7.4777*(pow(sol[5],3)+pow(sol[6],3))+0.7854*(sol[3]*pow(sol[5],2)+sol[4]*pow(sol[6],2)))+top;
    case 15: //10-bar-truss


	/*	truss(area, stress, displX, displY);         call truss FE to get lengths,
						 stresses and displacements

	init */
	for (i = 0; i < 12; i++)
		for (j = 0; j < 12; j++)
			a[i][j] = 0.;

	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
		       stiffmatrix[i][j] = 0.;

	area = (double*) malloc (dim * sizeof(double));
	for (i=0; i < dim; i++)
	{
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
	fX[0] = 0;
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
	for (m=0; m<neq;m++)
	{
		b[m] = 0.;
		for(k=0;k<neq;k++) a[m][k]=0.;
	}
	/*  >>>>>>> loop over all elements  <<<<<<< */
	for(m=0;m<num_elem;m++)
	{
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
		for (k = 0; k < 4; k++)
		{
			for(g=0; g<4; g++)
			{
				i = lm [k];
				j = lm [g];
				a[i][j] = a[i][j] + stiffmatrix[k][g];
			}
		}

	}
	/*end stiff*/

	/*displ*/
	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
		b[nn] = fY[n];
		b[nn-1] = fX[n];
	}

	for(n = 0; n < num_node; n++)
	{
		nn = (2 * n) + 1;
		switch(kode[n])
		{
			case 1:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-(a[i][nn-1]*fX[n]);
					a[i][nn-1]=0.;
					a[nn-1][i]=0.;
				}
				a[nn-1][nn-1]=1.00;
				b[nn-1]=fX[n];
				break;
			case 2:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-(a[i][nn]*fY[n]);
					a[i][nn]=0.;
					a[nn][i]=0.;
				}
				a[nn][nn]=1.00;
				b[nn]=fY[n];
				break;
			case 3:
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-a[i][nn-1]*fX[n];
					a[i][nn-1]=0.;
					a[nn-1][i]=0.;
				}
				a[nn-1][nn-1]=1.00;
				b[nn-1]=fX[n];
				for(i=0;i<neq;i++)
				{
					b[i]=b[i]-a[i][nn]*fY[n];
					a[i][nn]=0.;
					a[nn][i]=0.;
				}
				a[nn][nn]=1.00;
				b[nn]=fY[n];
				break;
		}
	};


	/*	symsol();*/
	neq1 = nl = 0;
	neq1=neq-1;
	nl = neq - 1;

	for (n=0;n<nl;n++)
	{
		if (a[n][n]<=0.)
		{
			printf("zero or negative main-diagonal %d\n",n);
			exit(2);
		}

		n1 = n +1 ;
		for (j=n1;j<neq;j++) {
			/*printf("l279: %f %f\n",a[n][j],a[n][n]);*/
			a[n][j]=a[n][j]/a[n][n];
		}
		for(i=n1;i<neq;i++)
		{
			if(a[n][i]==0.0)
				b[i]=b[i]-a[n][i]*b[n];
			else
			{
				for(j=i;j<neq;j++)
				{
					/*printf("l288: %f %f %f\n",a[i][j],a[i][n],a[n][j]);*/
					a[i][j]=a[i][j]-a[i][n]*a[n][j];
					a[j][i]=a[i][j];
				}
			}
			b[i]=b[i]-a[n][i]*b[n];
		}
		b[n]=b[n]/a[n][n];
	};

	m=neq1;
	b[m]=b[m]/a[m][m];
	for(n=0;n<nl;n++)
	{
		m1=m;
		m=m-1;
		for(j=m1;j<neq;j++) {
			b[m]=b[m]-b[j]*a[m][j];
			/*printf("l306: %f %f %f\n",b[m],b[j],a[m][j]);*/
		}
	}
	/*end symbol*/

	for(i=0;i<num_node;i++)
	{
		displX[i] = b[2*i];
		displY[i] = b[2*i+1];
	}

	/*printf("\n NODE    VEL-X     VEL-Y\n");
	for (m=0;m<num_node;m++)
		printf("%3d  %10.7f %10.7f\n", (m+1),displX[m],displY[m]);*/
	/*end displ*/


	/*stress*/
	/*printf("\nELEM I-NOD J-NOD    MEM-FOR     STRESS\n");*/
	for(m=0;m<num_elem;m++)
	{
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

		/*printf("%d    %d    %d  %10.4f  %10.4f\n", (m+1), (i+1), (j+1), p,strss[m]);*/
	}
	/*end stress*/

	for (i = 0; i < 10; i++) {
		/*weight of the truss*/
		somF += ro*area[i]*length[i];
	}
	for(i = 0; i < 10; i++) {
		float aux = fabs(stress[i]);
		somg1 += (aux > 25) ? aux - 25. : 0;          /*stress constraint (25 ksi)*/
	//		if (print_fitness == 1)
	//			printf("G1[%d] = %f\n", i, aux);
	}
	for(i = 0; i < 6; i++) {
		/*displacement constraint (2 in.)*/
		float aux = fabs(displX[i])/* * 1000.*/;
		somg2 += (aux > 2) ? aux - 2. : 0;
	//		if (print_fitness == 1)
	//			printf("G2X[%d] = %f\n", i, aux);
	}
	for(i = 0; i < 6; i++) {
		/*displacement constraint (2 in.)*/
		float aux = fabs(displY[i])/* * 1000.*/;
		somg2 += (aux > 2) ? aux - 2. : 0;
	//		if (print_fitness == 1)
	//			printf("G2Y[%d] = %f\n", i, aux);
	}
	//	if (print_fitness == 1) {
	//		printf("F %f G1 %f G2 %f\n", somF, somg1, somg2);
	//	}

        //func_eval[sp_index] += 1;

	return (somF+P*(somg1+somg2));                /*Transformed Function*/
	/////END 10-bar-truss
	default:
        printf("Info: Invalid function..\n") ;
        exit(0);
    }
}

void initPop(int func, int pop_size, int dim, double *best, double **pop, double *lb, double *ub, double *fo){
	int j,k;
	if(func!=10 && func!=11 && func!=12 && func!=13 && func!=14){
		for (j=0;j<pop_size;j++){//each individual
			fo[j]  = 0.0;
			for (k=0; k<dim;k++){ //each dimension of the individual
				best[k] = 0.0;
				pop[j][k] = randon(lb[0],ub[0]);
			}
		}
	}else{
		for (j=0;j<pop_size;j++){//each individual
			fo[j]  = 0.0;
			for (k=0; k<dim;k++){ //each dimension of the individual
				best[k] = 0.0;
				pop[j][k] = randon(lb[k],ub[k]);
			}
		}
	}
	
}

void mutualism_phase(double *y, int index_i, int dim, double **pop, double *best, double *ub, double *lb, double *fo, int func, int pop_size){
	int i;
	int bf1=0,bf2=0;	
	double *mutual;
	double *new_x_i,*new_x_j;
	double newfo_i,newfo_j;
	double *array_rand;
	int index_j=index_i;

	while(index_j==index_i){
		index_j=(int)randon(0,pop_size);
	}//pick xj	

	array_rand=(double*)malloc(dim*sizeof(double));
	new_x_i = (double*)malloc(dim*sizeof(double));
	new_x_j = (double*)malloc(dim*sizeof(double));	
	mutual = (double*)malloc(dim*sizeof(double));
	
	bf1=(int)randon(1,3);//benefit factor1
	bf2=(int)randon(1,3);//benefit factor2
	for(i=0;i<dim;i++){
		array_rand[i]=randon(0,1);
	}
	//cáculo Eqs(1) (2) (3)
	for(i=0;i<dim;i++){
		mutual[i]=(pop[index_i][i]+pop[index_j][i])/2;
		new_x_i[i]=pop[index_i][i]+(array_rand[i]*(best[i]-(mutual[i]*bf1)));
		new_x_j[i]=pop[index_j][i]+(array_rand[i]*(best[i]-(mutual[i]*bf2)));
		if(func!=10 && func!=11 && func!=12 && func!=13 && func!=14){
			if(new_x_i[i]>ub[0])new_x_i[i]=ub[0];
			if(new_x_j[i]>ub[0])new_x_j[i]=ub[0];
			if(new_x_i[i]<lb[0])new_x_i[i]=lb[0];
			if(new_x_j[i]<lb[0])new_x_j[i]=lb[0];
		}else{
			if(new_x_i[i]>ub[i])new_x_i[i]=ub[i];
			if(new_x_j[i]>ub[i])new_x_j[i]=ub[i];
			if(new_x_i[i]<lb[i])new_x_i[i]=lb[i];
			if(new_x_j[i]<lb[i])new_x_j[i]=lb[i];
		}
	}
	//
	newfo_i=objfunc(new_x_i, y, func, dim);//calculates fitness for new_x_i
	newfo_j=objfunc(new_x_i, y, func, dim);//calculates fitness for new_x_j

	if(fo[index_i]>=newfo_i){//greedy selection for xi
		for(i=0;i<dim;i++){
			pop[index_i][i]=new_x_i[i];
		}
		fo[index_i]=newfo_i;
	}

	if(fo[index_j]>=newfo_j){//greedy selection for xj
		for(i=0;i<dim;i++){
			pop[index_j][i]=new_x_j[i];
		}
		fo[index_j]=newfo_j;
	}

	free(array_rand);
	free(new_x_i);	
	free(new_x_j);
	free(mutual);		
}

void commensalism_phase(double *y, int index_i, int dim, double **pop, double *best, double *ub, double *lb, double *fo, int func, int pop_size){
	int i;
	double *array_rand;
	int index_j=index_i;
	double *new_x_i;
	double newfo_i;

	while(index_j==index_i){
		index_j=(int)randon(0, pop_size);
	}//pick xj	

	//alloc arrays
	array_rand=(double*)malloc(dim*sizeof(double));
	new_x_i = (double*)malloc(dim*sizeof(double));
	//

	//puts values for array_rand and new_x_i
	for(i=0;i<dim;i++){
		array_rand[i]=randon(-1,1);
		new_x_i[i]=pop[index_i][i]+(array_rand[i]*(best[i]-pop[index_j][i]));
		if(func!=10 && func!=11 && func!=12 && func!=13 && func!=14){
			if(new_x_i[i]<lb[0])new_x_i[i]=lb[0];
			if(new_x_i[i]>ub[0])new_x_i[i]=ub[0];
		}else{
			if(new_x_i[i]<lb[i])new_x_i[i]=lb[i];
			if(new_x_i[i]>ub[i])new_x_i[i]=ub[i];
		}
	}

	newfo_i=objfunc(new_x_i, y, func, dim);//calculates fitness for new_x_i

	if(fo[index_i]>=newfo_i){////greedy selection for xi
		for(i=0;i<dim;i++){
			pop[index_i][i]=new_x_i[i];
		}
		fo[index_i]=newfo_i;
	}
	free(array_rand);
	free(new_x_i);	
}

void parasitism_phase(double *y, int index_i, int dim, double **pop, double *ub, double *lb, double *fo, int func, int pop_size){
	double *parasite;
	double parasite_fo;
	int pick;
	int index_j=index_i;
	int i;

	parasite=(double*)malloc(dim*sizeof(double));//alloc parasite array	
	pick=(int)randon(0,dim);//chooses the dimension for change
	
	for(i=0;i<dim;i++){//copies xi for parasite
		parasite[i]=pop[index_i][i];
	}

	while(index_j==index_i){
		index_j=(int)randon(0,pop_size);
	}
	if(func!=10 && func!=11 && func!=12 && func!=13 && func!=14){
		parasite[pick]=randon(lb[0],ub[0]);//change the value of the dimension choosen
	}else{
		parasite[pick]=randon(lb[pick],ub[pick]);//change the value of the dimension choosen
	}
	parasite_fo=objfunc(parasite, y, func, dim);//calculates fitness for parasite
	
	if(fo[index_j]>=parasite_fo){//greedy selection between xj and parasite
		for(i=0;i<dim;i++){
			pop[index_j][i]=parasite[i];
		}
		fo[index_j]=parasite_fo;
	}
	free(parasite);	
}

void AvgStdDev(double *Avg,double *StdDev,double Var[], int run){
	int i;

	*Avg = 0;
	*StdDev = 0;
	for (i=0;i<run;i++){
		*Avg += Var[i];
	}
	*Avg /= run;

	for (i=0;i<run;i++){
		*StdDev += pow((*Avg-Var[i]),2);
	}
	*StdDev /= run;
   *StdDev = sqrt(*StdDev);
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
		fprintf(shellComands,"	plot \"dadosplot//dadosplot%s.txt\" using 1:3 title \"Media Fo\" w lines , \"dadosplot//dadosplot%s.txt\" using 1:2 title \"Best_Fo\" w lines",run,run);
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