
#ifndef LINA_STUFF
#define LINA_STUFF

#include <stdio.h>
#include <math.h>
#include <stdlib.h> 
#include <string.h>

///////////////////////////////////////// public structs ////////////////////////////////////////////////////////////////////////

typedef struct CSR{
	int rows;
	int* row_start;
	int* indices;
	double* elements;
}CSR_matrix;

typedef struct CSC{
	int cols;
	int* col_start;
	int* indices;
	double* elements;
}CSC_matrix;

//////////////////////////////////////// public functions ///////////////////////////////////////////////////////////////////////

// print functions
void print_vector(double* V,int n,const char* Name);
void print_list(int* V,int n,const char* Name);
void print_CSR(CSR_matrix* A,const char* Name);
void print_CSC(CSC_matrix* A,const char* Name);
void print_CSR_matrix_market(CSR_matrix* A,const char* Name);
void print_CSC_matrix_market(CSC_matrix* A,const char* Name);

// vector functions
double* zero_vector(int n);
void clear_vector(double* V,int n);

// matrix functions
void free_CSR_matrix(CSR_matrix** const A);
void free_CSC_matrix(CSC_matrix** const A);
void CSR_scale_matrix(CSR_matrix* A,double factor);
CSR_matrix* matrix_product_to_CSR(CSR_matrix* A,CSC_matrix* B);


// miscellaneous functions
int* zero_int_list(int n);

#endif
