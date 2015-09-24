#ifndef SOSFUNCTIONS_H
#define SOSFUNCTIONS_H

typedef struct{
	int tid;
	int inicio;
	int fim;
	//double *best_th;
	//double bestfo_th;
	//int best_index_th;
}slice;


/* Control Parameters of the search algorithm*/
int POP_SIZE;  /* The number of candidate solutions*/
int MAX_ITER; /*The number of iterations*/
int FUNCTION;
int CORES;
//Problem definitions
int DIM;//number of problem variables
double *lb;//lower bound of the variables
double *ub;//upper bound of the variables
int RUN;/*Algorithm can run many times in order to see its robustness*/
double bestfoRUN;
//Global variables
double **pop;//[POP_SIZE][DIM];  //population of candidate solutions.
double *fo;//[POP_SIZE];      //objective function value.
double *best;//[DIM];           //best solution found
double bestfo;//best fo value
int best_index;//index for the best solution
double *y; //constraint functions
double *c; //penallity function constant
int rest;
int num_fit_eval;
int num_iter;

slice *argThread;



//Functions declarations
double randon(double inferior, double superior);
void freeArrays();
int GetParameters(char **argv);
void showParameters();
void prepararObjFunc();
double constr(double *sol);
double objfunc(double *sol, int cond);
void *th_init_pop(void *arg);
void initPop();
void mutualism_phase(int index_i, double **pop_, double *best_, double *fo_);
void commensalism_phase(int index_i, double **pop_, double *best_, double *fo);
void parasitism_phase(int index_i, double **pop_, double *fo);
int AvgStdDev(double *Avg,double *StdDev,double Var[]);
void plot(FILE *shellComands, char *run);
char *converteDecChar(char *strf, int dec);
int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer);
void sos_iter();
void alloc_slice_pointers();
void cpy_best();
void cpy_slice_pointers();
void join_pop();
void join_best();
void free_slice();

#endif
