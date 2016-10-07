#include "easyPCG.h"
#include "linalg_stuff.h"
#include <stdlib.h>
#include <time.h>

/* 
 * This is a test program for the PCG algorithm  
 * It computes the solution for a Poisson equation in 2D on a square domain
 * with a point source in the middle. The output is a textfile with the 
 * iteration number and the computation time (using C's clock()-function) for 
 * 10 different system sizes. This is done for all available preconditioners 
 * except for the user-defined. 
 * Copy "easyPCG.c, easyPCG.h, linalg_stuff.c, linalg_stuff.h" in the folder where the testPCG 
 * program is compiled. Then run make. 
 * 
 * 
 * The arguments passed to the program are the following:
 * 
 * #1: system size (i.e. node number)
 * #2: absolute name to output file
 * #3: absolute filename of the solution vector if desired (optionally)
 * 
 * Examples are:
 * 	./testPCG 5000 /Home/<MyDir>/PCGstat /Home/MyDir/solution
 *  or 
 * ./testPCG 5000 /Home/<MyDir>/PCGstat 
 * without solution
 * 
 * In the first example the programm will output the solutions 'solution_1,solution_2, ...'
 * in /Home/<MyDir>
 * 
 * 
 * The structure of the PCGstat output file is:
 * 
 * col1 :system size ; col2: iterations with precond #1 ; col3: execution time with precond #1 ; col4; iterations with precond #2; etc.
 * 
 * 
 */


CSR_matrix* Laplace1D_3point(int size){
	int i;
	
	CSR_matrix* Res = (CSR_matrix*)malloc(sizeof(CSR_matrix));
	Res->rows = size;
	Res->row_start = (int*)malloc((size+1)*sizeof(int));
	Res->indices = (int*)malloc((3*(size-2)+4)*sizeof(int));
	Res->elements = (double*)malloc((3*(size-2)+4)*sizeof(double));
	
	int pos = 0;
	for (i=0;i<size;i++){
		Res->row_start[i] = pos;
		
		if (i>0){
			Res->indices[pos] = i-1;
			Res->elements[pos++] = 1.0;
		}
		Res->indices[pos] = i;
		Res->elements[pos++] = -2.;
		if (i<size-1){
			Res->indices[pos] = i+1;
			Res->elements[pos++] = 1.0;		
		}
	}
	Res->row_start[i] = pos;
	
	return Res;
}

CSR_matrix* Laplace2D_5point(int* size){
	int i,j;
	
	int width = (int)round(sqrt(*size));
	*size = width*width;
	
	CSR_matrix* Res = (CSR_matrix*)malloc(sizeof(CSR_matrix));
	Res->rows = *size;
	Res->row_start = (int*)malloc((*size+1)*sizeof(int));
	Res->indices = (int*)malloc(5*(*size)*sizeof(int));
	Res->elements = (double*)malloc(5*(*size)*sizeof(double));
	
	int pos = 0;
	int var = 0;
	for (i=0;i<width;i++){
		for (j=0;j<width;j++){
			Res->row_start[var] = pos;
			
			if (i>0){
				Res->indices[pos] = var-width;
				Res->elements[pos++] = 1.;
			}
			
			if (j>0){	
				Res->indices[pos] = var-1;
				Res->elements[pos++] = 1.;
			}
			
			Res->indices[pos] = var;
			Res->elements[pos++] = -4.;
			
			if (j<width-1){	
				Res->indices[pos] = var+1;
				Res->elements[pos++] = 1.;
			}				
			
			if (i<width-1){
				Res->indices[pos] = var+width;
				Res->elements[pos++] = 1.;
			}
			
			var++;
		}
	}
	Res->row_start[var] = pos;
	
	Res->indices = (int*)realloc(Res->indices,pos*sizeof(int));
	Res->elements = (double*)realloc(Res->elements,pos*sizeof(double));

	return Res;
}

double* SourceTerm1D(int size,double left_val,double right_val){
	int i;
	
	double* Res = (double*)calloc(size,sizeof(double));
	
	int n1 = size/4;
	int n2 = 3*n1;
	
	Res[n1] = 1.;
	Res[n2] = -5.;	
	
	for (i=0;i<size;i++) Res[i] /= (double)(size-1)*(size-1);
	
	Res[0] = -left_val;
	Res[size-1] = -right_val;
	
	return Res;
}

double* SourceTerm2D(int size,double bound_val){
	int i,j;
	double x,y,f;
	
	int width = (int)sqrt(size);
	double* Res = (double*)calloc(size,sizeof(double));
	
	for (i=0;i<width;i++){
		x = (double)i/(width-1);
		for (j=0;j<width;j++){
			y = (double)j/(width-1);
			f = exp(-20.*((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)));
			Res[i*width+j] = (double)f/((width-1)*(width-1));			
		}
	}
		
	
	for (i=0;i<width;i++){
		Res[i] += bound_val;
		Res[i*width] += bound_val;
		Res[(i+1)*width-1] += bound_val;
		Res[width*(width-1)+i] += bound_val;		
	}
	
	return Res;
}

int main(int argc, char* argv[]){
	
	const int number_of_sols = 10;												// number of simulations with different system sizes
	const int max_iter = 1000;													// maximum number of allowed iterations 
	const double tol = 1e-11;													// required residuum reduction
	
	int i,max_size,mode,iter,size;
	char line[512],strval[512];
	
	clock_t t,t0; 
	
	char* solution_out = NULL;
	FILE* file = NULL;
	
	if (argc==1){
		max_size = 1000;
		printf("no arguments given -> use some standard standard value: system size = %d\n",max_size);
		printf("no output wil be generated\n\n");
	}
	else{
		if (argc<3){
		 printf("error in %s: not enough arguments given\n",__FILE__);
		 exit(EXIT_FAILURE);
		}
	
		max_size = atoi(argv[1]);
		if (max_size<=0){
			printf("error in %s: invalid argument #1, found: %d, required: positive integer\n",__FILE__,max_size);
			exit(EXIT_FAILURE);		
		}
	
		file = fopen(argv[2],"w");
		if (file==NULL){
			printf("error in %s: could not open file %s\n",__FILE__,argv[2]);
			exit(EXIT_FAILURE);
		}	
	
		if (argc>3){
			solution_out = (char*)malloc(512*sizeof(char));
			strcpy(solution_out,argv[3]);
		}
	}
	
	printf("\nstart\n");
	fflush(stdout);
	for (i=1;i<=number_of_sols;i++){
		
		size = (max_size / number_of_sols)*i;
	
		/*
		// Creates a Poisson-equation system in 1D with Dirichlet boundary conditions X(left)=1, X(right)=1
		CSR_matrix* A = Laplace1D_3point(max_size);
		double* B = SourceTerm1D(max_size,1.,1.);
		*/
	
		// Creates a Poisson-equation system in 2D with Dirichlet boundary conditions X(bound)=0 and a point source in the center
		CSR_matrix* A = Laplace2D_5point(&size);
		double* B = SourceTerm2D(size,0.);	
		CSR_scale_matrix(A,-1.);
		double* Sol = zero_vector(size);
		
		printf("solving with %d nodes ...\n",size);
		fflush(stdout);
		sprintf(line,"%d",size);
		for (mode=PCG_NONE;mode<=PCG_ICHOL;mode++){
			if ((preconditionerID)mode!=PCG_USER_DEFINED){								
				PCG_set_preconditioner_mode(A,(preconditionerID)mode);							
				t0 = clock();	
				clear_vector(Sol,size);
				iter = PCG_solve(Sol,B,size,tol,max_iter);
				t = clock();
				if (iter>=max_iter){
					printf("PCG with preconditioner no. %d did not converge %d iterations with tolerance %e\n",mode,max_iter,tol);
				}
				sprintf(strval,"\t%d\t%f",iter,(double)(t-t0)/CLOCKS_PER_SEC);
				strcat(line,strval);												
			}
		}
		if (file!=NULL) fprintf(file,"%s\n",line);
		PCG_clean();
		free_CSR_matrix(&A);
		free(B);
		
		if (solution_out!=NULL){
			sprintf(line,"%s_%d",solution_out,i);
			print_vector(Sol,size,line);
		}
		free(Sol);
	}
	
	printf("\nfinished -> test successfull !\n");
	if (file!=NULL) fclose(file);
	if (solution_out!=NULL) free(solution_out);
	return 0;
}
