#ifndef SOSFUNCTIONS_H
#define SOSFUNCTIONS_H

//Functions declarations

	double randon(double inferior, double superior);
	//void AllocArrays(int pop_size, int dim, double **pop, double *fo, double *best, double *ub, double *lb);
	void freeArrays(int pop_size, double **pop, double *fo, double *best, double *ub, double *lb);
	int GetParameters(char **argv, int *run, int *max_iter, int *pop_size, int *dim, int *function);
	void showParameters(int func, int run, int max_iter, int pop_size, int dim);
	void prepararObjFunc(int func, double *lb, double *ub);
	double constr(double sol[], int func, double *y);
	double objfunc(double sol[], double y[], int func, int dim);
	void initPop(int func, int pop_size, int dim, double *best, double **pop, double *lb, double *ub, double *fo);
	void mutualism_phase(double *y, int index_i, int dim, double **pop, double *best, double *ub, double *lb, double *fo, int func, int pop_size);
	void commensalism_phase(double *y, int index_i, int dim, double **pop, double *best, double *ub, double *lb, double *fo, int func, int pop_size);
	void parasitism_phase(double *y, int index_i, int dim, double **pop, double *ub, double *lb, double *fo, int func, int pop_size);
	void AvgStdDev(double *Avg,double *StdDev,double Var[], int run);
	void plot(FILE *shellComands, char *run);
	char *converteDecChar(char *strf, int dec);
	int ffscanf(char *fieldname, FILE *fp, char *format, void *inbuffer);

#endif