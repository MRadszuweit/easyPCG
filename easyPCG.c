
#include "easyPCG.h"

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////* local structs *////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Matrix type in Compressesed Sparse Row (CSR) format
 * 
 * see https://en.wikipedia.org/wiki/Sparse_matrix for details of the data format
 * 
 */
typedef struct CSR{
	int rows;						/**< number of rows */ 
	int* row_start;					/**< indices that point to the beginning of each row in indices */
	int* indices;					/**< list of all indices */
	double* elements;				/**< list of all values */
}CSR_matrix;

/**
 * Matrix type in Compressesed Sparse Column (CSC) format
 * 
 * see https://en.wikipedia.org/wiki/Sparse_matrix for details of the data format
 * 
 */
typedef struct CSC{
	int cols;						/**< number of columns */ 
	int* col_start;					/**< indices that point to the beginning of each row in indices */
	int* indices;					/**< list of all indices */
	double* elements;				/**< list of all values */
}CSC_matrix;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////* local (global) variables *//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

static preconditionerID precon_mode = PCG_NONE;												// Preconditioner mode ID
static void (*precon)(double* R,double* X,int n) = NULL;								// approximates (A^-1).X and stores in R
static void (*user_def_mult)(double* R,double* X,int n) = NULL;							// user-defined multiplcation routine R=A.X
static CSR_matrix* default_A = NULL;													// standard matrix in Compressed Row Storage (CSR) format
static CSC_matrix* IL = NULL;															// if ILU is used: incomplete Lower Triangular matrix
static CSC_matrix* ILT = NULL;															// if ILU is used: the transposed of the above

static double pcg_residuum;																// the final residuum is stored in this variable


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////* local helper functions *//////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//// list functions
/**
 * \brief Allocates a list of integers initialized by zero
 * \param int size: size of list
 * \return int* : pointer to allocated list
 */
static int* zero_int_list(int n){														
	return (int*)calloc(n,sizeof(int));
}





/**
 * \brief Inserts an integer into an ordered list
 * \param int* list: pointer to list
 * \param double* vals: pointer to list of double values associated with each integer
 * \param int size: number of entries in both lists
 * \param int index: index that is inserted 
 * \param double value: double value that will be inserted in the double-valued list
 * 
 * This method operates on a row/column of a matrix with ordered indices. 
 * The row/column is reporesented by a list of column/row (list) indices and the 
 * associated matrix values (vals). Both, the index and the element values are inserted, 
 * such that the ordering is preserved (in ascending order).
 * The method is not optimizeed yet but is sufficient for small row/column sizes (~10)
 */
static void insert_int_list(int* list,double* vals,int size,int index,double value){	// inserts element in ordered integer list 
	int i = 0;
	while (i<size){
		if (index<list[i]) break;
		i++;
	}
	int d = size-i;
	if (d>0){
		memmove(&(list[i+1]),&(list[i]),d*sizeof(int));
		memmove(&(vals[i+1]),&(vals[i]),d*sizeof(double));
	}
	list[i] = index;	
	vals[i] = value;
}


//// vector related functions
/**
 * \brief Allocates a list of doubles initialized by zero
 * \param int size: dimension of vector
 * \return double* : pointer to allocated list
 */
static double* zero_vector(int size){														
	return (double*)calloc(size,sizeof(double));
}




/**
 * \brief Returns a copy of a vector 
 * \param double* X: pointer to the source vector
 * \param int n: dimension of vector
 * \return double* : pointer to the newly allocated copy of vector X
 */
static double* copy_vector(double* X,int n){						
	double* R  = (double*)malloc(n*sizeof(double));
	memmove((void*)R,(void*)X,n*sizeof(double));
	return R;
}




/**
 * \brief Scales a vector X <- factor*X
 * \param double* X: pointer to the vector te be scaled
 * \param double factor: scaling factor
 * \param int size: dimension of vector
 */
static void vector_scale(double* X,double factor,int size){								
	int i;
	for (i=0;i<size;i++){
		X[i] *= factor;
	}
}




/**
 * \brief Adds a scaled vector X <- X+factor*Y
 * \param double* X: pointer to the vector to be manipulated
 * \param double* Y: pointer to the vector that will be added
 * \param double factor: scaling factor
 * \param int size: dimension of vector
 */
static void vector_addmult(double* X,double* Y,double factor,int size){	
	int i;
	for (i=0;i<size;i++){
		X[i] += factor*Y[i];
	}
}




/**
 * \brief Computes the euklidian scalar product of two vectors
 * \param double* X: first argument
 * \param double* Y: second argument
 * \return double: scalar product <X,Y>
 */
static double vector_scalar(double* X,double* Y,int size){							
	int i;
	double sum = 0;
	for (i=0;i<size;i++){
		sum += X[i]*Y[i];
	}
	return sum;
}




/**
 * \brief Computes the euklidian norm of a vector
 * \param double* X: pointer to vector 
 * \return double: euklidean norm ||x|| = sqrt(<X,X>)
 */
static double euklid_norm(double* X,int size){	
	return sqrt(vector_scalar(X,X,size));
}




//// matrix related functions
/**
 * \brief frees all memory of a CSR-matrix struct
 * \param CSR_matrix** A: pointer to the pointer where the matrix is stored 
 */
static void free_CSR_matrix(CSR_matrix** A){	
	if (*A!=NULL){
		if ((*A)->row_start!=NULL) free((*A)->row_start);
		if ((*A)->indices!=NULL) free((*A)->indices);
		if ((*A)->elements!=NULL) free((*A)->elements);
		free(*A);
		*A = NULL;
	}
}




/**
 * \brief Frees all memory of a CSC-matrix struct
 * \param CSC_matrix** A: pointer to the pointer where the matrix is stored 
 */
static void free_CSC_matrix(CSC_matrix** A){			
	if (*A!=NULL){
		if ((*A)->col_start!=NULL) free((*A)->col_start);
		if ((*A)->indices!=NULL) free((*A)->indices);
		if ((*A)->elements!=NULL) free((*A)->elements);
		free(*A);
		*A = NULL;		
	}
}




/**
 * \brief Computes the tranpose of a CSR-matrix as a CSC-matrix
 * \param CSR_matrix* A: pointer the the CSR-matrix to be tranposed
 * \return CSC_matrix*: pointer to the newly allocated CSC-matrix that is the transposed of A (very fast)
 */
static CSC_matrix* transpose_CSR_to_CSC(CSR_matrix* A){								
	int n = A->rows;											
	int N = A->row_start[n];
	CSC_matrix* Res = (CSC_matrix*)malloc(sizeof(CSC_matrix));
	Res->cols = n;
	Res->col_start = zero_int_list(n+1);
	Res->indices = zero_int_list(N);
	Res->elements = zero_vector(N);
	
	memmove((void*)Res->col_start,(void*)A->row_start,(n+1)*sizeof(int));
	memmove((void*)Res->indices,(void*)A->indices,N*sizeof(int));
	memmove((void*)Res->elements,(void*)A->elements,N*sizeof(double));	
	return Res;
}




/**
 * \brief Solves the a linear equation A.X=B, when the LU-decomposition of A=L.U is given
 * \param CSC_matrix* L: pointer to the lower triangular matrix of the LU decompostion
 * \param CSC_matrix* U: pointer to the upper triangular matrix of the LU decompostion
 * \param double* X: pointer to vector where the result is stored in
 * \param double* B: pointer to the right side of the equation A.X=B, A=L.U
 * 
 * This method solves the the equation A.X=B, when the LU decompostion A=L.U is known 
 * via forward/backward substitution. Both matrices must be given in column-format (CSC), 
 * in order to better parallelize the task later. 
 */
static void CSC_triangular_LU_invert(CSC_matrix* L,CSC_matrix* U,double* X,double* B){
	int j,I,J,start,end;
	double x;
	
	int n = L->cols;
	X = memmove((void*)X,(void*)B,n*sizeof(double));
	
	// L.Y=B
	for (I=0;I<n;I++){
		start = L->col_start[I];
		x = X[I]/L->elements[start++];
		X[I] = x;
		for (j=start;j<L->col_start[I+1];j++){
			J = L->indices[j];
			X[J] -= L->elements[j]*x;
		}
	}	
	
	// U.X=Y	
	for (I=n-1;I>=0;I--){
		end = U->col_start[I+1]-1;
		x = X[I]/U->elements[end--];
		X[I] = x;
		for (j=end;j>=U->col_start[I];j--){
			J = U->indices[j];
			X[J] -= U->elements[j]*x;
		}
	}	
}




/**
 * \brief Converts a CSR-matrix to the CSC storage format
 * \param CSR_matrix* A: pointer the the matrix in row-storage format (CSR) the will be converted
 * \param int cols: column number of the matrix A
 * \return CSC_matrix*: copy of the CSR matrix in column-storage format (CSC) 
 * 
 * Since the column number is not explicitly stored in the CSR format, it must be passed
 * to the method. For invertable matrices it is usually #rows = #colums, which makes it simple. 
 */
static CSC_matrix* CSR_to_CSC_matrix(CSR_matrix* A,int cols){
	
	int i,j,col;
	
	int n = A->rows;
	int N = A->row_start[n];
	int mean_bandwidth = N/n+1;	
	int* col_ptr = zero_int_list(cols);
	int* mem_sizes = zero_int_list(cols);
	int** Col_indices = (int**)malloc(cols*sizeof(int*));
	double** Col_elements = (double**)malloc(cols*sizeof(double*));
	for (i=0;i<cols;i++){
		mem_sizes[i] = mean_bandwidth;
		Col_indices[i] = zero_int_list(mean_bandwidth);
		Col_elements[i] = zero_vector(mean_bandwidth);
	}
	
	// construct (ordered) column lists
	for (i=0;i<n;i++){
		for (j=A->row_start[i];j<A->row_start[i+1];j++){
			col = A->indices[j];
			if (col>=cols){
				printf("error in function %s in file %s: index %d exceeds given column number %d -> abort\n",__func__,__FILE__,col,cols);
				exit(EXIT_FAILURE);
			}						
			// use index-ordered insert here 
			insert_int_list(Col_indices[col],Col_elements[col],col_ptr[col],i,A->elements[j]);						
			col_ptr[col]++;
			if (col_ptr[col]==mem_sizes[col]-1){
				mem_sizes[col] += mean_bandwidth;
				Col_indices[col] = (int*)realloc(Col_indices[col],mem_sizes[col]*sizeof(int));
				Col_elements[col] = (double*)realloc(Col_elements[col],mem_sizes[col]*sizeof(double));
			}
			
		}
		
	}
	
	CSC_matrix* Res = (CSC_matrix*)malloc(sizeof(CSC_matrix));
	Res->cols = cols;
	Res->col_start = zero_int_list(cols+1);
	Res->indices = zero_int_list(N);
	Res->elements = zero_vector(N);
	
	// assumes ordered columns
	int Res_ptr = 0;
	for (i=0;i<cols;i++){		
		for (j=0;j<col_ptr[i];j++){
			Res->indices[Res_ptr] = Col_indices[i][j];
			Res->elements[Res_ptr] = Col_elements[i][j];
			Res_ptr++;
		}		
		Res->col_start[i+1] = Res_ptr;
	}
	
	// clean
	for (i=0;i<cols;i++){
		free(Col_indices[i]);
		free(Col_elements[i]);
	}
	free(Col_indices);
	free(Col_elements);
	free(mem_sizes);
	free(col_ptr);
	
	return Res;
}




/**
 * \brief Helper function for Cholesky decomposition
 * \param CSR_matrix* L: Pointer where the lower triangular factor of the decomposition is stored
 * \param int row1: index of first row
 * \param int row2: index of second row
 * \param int max_ind: index limiter for the row product
 * \return double: row product inside matrix L: \sum_{k<max_ind} L_{row1,k}*L_{row2,k}
 */
static double CSR_row_product(CSR_matrix* L,int row1,int row2,int max_ind){	
	
	int start_i = L->row_start[row1];
	int end_i = L->row_start[row1+1];
	int start_j = L->row_start[row2];	
	int end_j = L->row_start[row2+1];
	
	int i = start_i;
	int j = start_j;
	int I = L->indices[i];
	int J = L->indices[j];
	
	double sum = 0;
	while(i<end_i){
		I = L->indices[i];
		if (I>=max_ind) break;				
		while(j<end_j){
			J = L->indices[j];
			if (J>=I) break;
			j++;
		}
		if (J==I) sum += L->elements[i]*L->elements[j];
		i++;
	}	
	return sum;
}




/**
 * \brief Gives a diagonal element of matrix A
 * \param CSR_matrix* A: pointer to matrix
 * \param int row: row index of the diagonal element 
 * \return double: diagonal element A_{row,row} or zero if row exceeds the row number of the matrix
 */
static double CSR_get_diag_element(CSR_matrix* A,int row){
	if (row>=A->rows) return 0;
	
	int j = A->row_start[row];
	int end = A->row_start[row+1];
	int J = A->indices[j];
	while(j<end){
		J = A->indices[j];
		if (J>=row) break;
		j++;
	}
	if (J==row) return A->elements[j]; else return 0;
}




/**
 * \brief Computes the incomplete Cholesky decomposition of a matrix
 * \param CSR_matrix* A: pointer to the matrix to be factorized
 * \param CSR_matrix* L: pointer where the resulting factor L.L^T=A is stored. Allocate in uninitialized CSR struct and pass the pointer.
 * 
 * This method computes the incomplete Cholesky decomposition A=L.L^T. The matrix to be factorized must 
 * be symmetric and positive definite to ensure existence of the factorization. For negative definete symmetric matrices 
 * multiply by -1 to obtain a positive definite matrix before decomposition. 
 */
void CSR_incompleteCholeski(CSR_matrix* A,CSR_matrix* L){
	int start,end,j,k,I,J;
	double Lii,Lij;
	
	int n = A->rows;
	int N = A->row_start[n];
	
	// create pattern of L (lower diagonal)
	L->rows = n;
	L->row_start = zero_int_list(n+1);
	L->indices = zero_int_list(N);
	L->elements = zero_vector(N);
	L->row_start[0] = 0;
	
	for (I=0;I<n;I++){
		start = A->row_start[I];
		end = A->row_start[I+1];							
		
		j = start;		
		k = L->row_start[I];
		while(j<end){	
			J = A->indices[j++];
			if (J>I) break;
			L->indices[k++] = J;
		}				
		L->row_start[I+1] = k;				
	}
	int M = L->row_start[n];
	L->indices = (int*)realloc(L->indices,M*sizeof(int));
	L->elements = (double*)realloc(L->elements,M*sizeof(double));
	
	
	int* col_ptrs_L = zero_int_list(n);
	memmove(col_ptrs_L,L->row_start,n*sizeof(int));
	
	for (I=0;I<n;I++){
		Lii = CSR_get_diag_element(A,I)-CSR_row_product(L,I,I,I);				
		if (Lii<=0){
			printf("error in function %s in file %s: improper values found at row %d-> abort\n Maybe you should check if your matrix is symmetric and positive-definite\n",__func__,__FILE__,I);
			exit(EXIT_FAILURE);
		}	
		Lii = sqrt(Lii);
		L->elements[col_ptrs_L[I]++] = Lii;
		
		start = A->row_start[I];
		j = A->row_start[I+1]-1;		
		while(j>=start){
			J = A->indices[j];
			if (J<=I) break;			
			Lij = (A->elements[j]-CSR_row_product(L,I,J,I))/Lii;
			L->elements[col_ptrs_L[J]++] = Lij;						
			j--;
		}		
	}
	
	free(col_ptrs_L);
}




/**
 * \brief Computes the matrix-vector product
 * \param double* R: pointer to vector where the product is stored. The memory must be allovcated before passing to the method.
 * \param CSR_matrix* A: pointer to the mapping matrix
 * \param double* X: pointer to vector that is multiplied by A
 * 
 * The memory where the result R is stored must be allocated before passed to this function. 
 * The dimension of R must be equal the the row number of the matrix A
 */
static void CSR_matrix_mult(double* R,CSR_matrix* A,double* X){		
		
	int* ptr_start = A->row_start;
	int* ptr_ind = A->indices;
	double* ptr_mat = A->elements;
	double* ptr_res = R;
	
	int i = 0;
	int n = A->row_start[A->rows];
	while (i<n){
		ptr_start++;
		*ptr_res = 0;						// if removed, then A.x is added to R
		while (i<*ptr_start){
			*ptr_res += (*ptr_mat)*X[*ptr_ind]; 
			ptr_ind++;
			ptr_mat++;
			i++;
		}
		ptr_res++;		
	}
}




/**
 * \brief Multiplication routine
 * \param double* R: pointer where the result of multiplicatione with the system matrix is stored
 * \param double* X: pointer to vector that is multiplicated
 * \param int n: dimension of the vectors
 */
static void PCG_mult(double* R,double* X,int n){
	if (user_def_mult!=NULL){
		user_def_mult(R,X,n);
	}
	else if (default_A!=NULL){		
		CSR_matrix_mult(R,default_A,X);
	}
	else{
		printf("error in function %s in file %s: no system matrix and no user-defined multiplication routine set -> abort",__func__,__FILE__);
		exit(EXIT_FAILURE);
	}
}




////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////7////////////////////* preconditioners */////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * \brief Dummy preconditioner
 * \param double* R: pointer to memory where the result ist stored
 * \param double* X: pointer to vector in which the preconditioner is applied
 * \param int n: dimension of vectors
 * 
 * This is a dummy preconditioner, i.e. no preconditioning is applied. 
 * The memory where R points to must be allocated before passing. 
 */   
void dummy_precon(double* R,double* X,int n){		
	memmove((void*)R,(void*)X,n*sizeof(double));
}




/**
 * \brief Simple Jacobi preconditioner
 * \param double* R: pointer to memory where the result ist stored
 * \param double* X: pointer to vector in which the preconditioner is applied
 * \param int n: dimension of vectors
 * 
 * This is the simplest preconditioner M.R=X. The preconditioning operation X=M^{-1}.R is performed 
 * by setting M=diagonal(A), where A is the system matrix. The system matrix must be stored in the 
 * global variable "default_A". The memory where R points to must be allocated before passing. 
 * 
 */
void CSR_Jacobi_precon(double* R,double* X,int n){					
	if (default_A!=NULL){
		int i,j;
		double aii;
		
		for (i=0;i<default_A->rows;i++){
			for (j=default_A->row_start[i];j<default_A->row_start[i+1];j++){
				if (default_A->indices[j]==i){
					aii= default_A->elements[j];
					if (aii!=0) R[i] = X[i]/aii;
				} 
			}
		}
	}
}




/**
 * \brief Gauss-Seidel preconditioner
 * \param double* R: pointer to memory where the result ist stored
 * \param double* X: pointer to vector in which the preconditioner is applied
 * \param int n: dimension of vectors
 * 
 * maybe useless
 */
void CSR_GausSeidel_precon(double* R,double* X,int n){	
	if (default_A!=NULL){
		int i,j;
		double sum;		
		for (i=0;i<default_A->rows;i++){
			sum = X[i];
			for (j=default_A->row_start[i];j<default_A->row_start[i+1];j++){
				if (default_A->indices[j]<i) sum -= default_A->elements[j]*R[default_A->indices[j]];
				if (default_A->indices[j]==i){
					R[i] = sum/default_A->elements[j];								
					break;
				}
			}
		}				
	}
}




/**
 * \brief Incomplete-Cholesky preconditioner
 * \param double* R: pointer to memory where the result ist stored
 * \param double* X: pointer to vector in which the preconditioner is applied
 * \param int n: dimension of vectors
 * 
 * This preconditioner assumes M=L.L^{T}, where L is the incomplete Cholesky factor of the 
 * system matrix. The preconditioning is basically forward/backward solving of a triangular system. 
 * The factor L and its transpose L^{T} must be stored in the global variables IL and ILT, respectively.
 * The memory where R points to must be allocated before passing. 
 */
void CSC_IC_precon(double* R,double* X,int n){								// incomplete Cholesky decomposition
	if (IL==NULL || ILT==NULL){
		printf("error in funciton %s in file %s: no ILU decomposition found -> abort\n",__func__,__FILE__);
		exit(EXIT_FAILURE);
	}	
	CSC_triangular_LU_invert(IL,ILT,R,X);
}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////* public functions *///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/** 
 * \brief Sets the system matrix A, for which A.X=B will be solved
 * \param void* A: pointer to the system matrix A
 * 
 * The system matrix is passed as a void pointer but make sure that the struct 
 * it points to has the same structure as CSR_matrix in external libraries. 
 */
void PCG_set_default_matrix(void* A){
	default_A = (CSR_matrix*)A;
}




/**
 * \brief Sets the preconditioner
 * \param void* System_matrix: pointer to the system matrix A, for which A.X=B will be solved
 * \param preconditionerID mode: preconditioner-mode ID
 * 
 * This method sets the preconditioner that will be used for the PCG algorithm. Since for some preconditioners 
 * precomputations must be performed, the system matrix can be passed, too. 
 * It is passed as a void pointer but make sure that the struct it points to, has the same structure as CSR_matrix 
 * in external libraries. 
 * 
 * If System_matrix is set to NULL, the default matrix does not change. This can make sense, if one wants 
 * to change the preconditioner method but not the system to be solved. 
 * 
 * If the System_matrix is not NULL the global variable default_A is set and it is not necessary to call 
 * "PCG_set_default_matrix" any more. 
 * 
 * There are currently 4 preconditioner modes avialable: 
 * 
 * 1. mode==PCG_NONE: no preconditioner is used
 * 2. 		PCG_USER_DEFINED: user defined preconditioner
 * 							If this mode is used, the method "PCG_set_precon" must be called to 
 * 							set the preconditioner.
 * 3. 		PCG_JACOBI: Jacobi preconditioner 
 * 4.		PCG_ICHOL: Incomplete Cholesky preconditioner 
 * 
 */
void PCG_set_preconditioner_mode(void* System_matrix,preconditionerID mode){	
		
	CSR_matrix* A = (CSR_matrix*)System_matrix;
	precon_mode = mode;
	if (A!=NULL) default_A = A;
	
	switch(precon_mode){
		case PCG_NONE:
			precon = dummy_precon;
			break;
		case PCG_USER_DEFINED:
			precon = NULL;
			break;
		case PCG_JACOBI:
			precon = CSR_Jacobi_precon;
			break;
		case PCG_ICHOL: 
			precon = CSC_IC_precon;
			if (A==NULL){
				if (default_A!=NULL) A = default_A;
				else{
					printf("warning: ILU preconditioner used but no system matrix defined -> skip\n");
					return;
				}
			}						
			CSR_matrix* Row_L = (CSR_matrix*)malloc(sizeof(CSR_matrix));			
			CSR_incompleteCholeski(default_A,Row_L);			
			
			if (IL!=NULL) free_CSC_matrix(&IL);
			if (ILT!=NULL) free_CSC_matrix(&ILT);
			IL = CSR_to_CSC_matrix(Row_L,default_A->rows);
			ILT = transpose_CSR_to_CSC(Row_L);
			
			free_CSR_matrix(&Row_L);			
			break;
		default:
			precon = dummy_precon;
			precon_mode = PCG_NONE;
			printf("warning: no valid preconditioner ID given -> NONE preconditioner set\n");
	}
}




/**
 * \brief Cleans internally allocated numerical objects used by the PCG algorithm
 * 
 */
void PCG_clean(){
	if (IL!=NULL) free_CSC_matrix(&IL);
	if (ILT!=NULL) free_CSC_matrix(&ILT);
}




/**
 * \brief Sets a user-defined preconditioner
 * \param void (*Precon)(double* R,double* X,int n): pointer to the preconditioner method.
 * 
 * This function must be called when the preconditioner mode "USER_DEFINED" is set. 
 * The Argument Precon must perform a preconditioner operation on X. The result M^{-1}.X must stored 
 * at pointer R. n denotes the dimension of the vectors. 
 */
void PCG_set_precon(void (*Precon)(double* R,double* X,int n)){
	precon = Precon;
}




/**
 * \brief Sets a user-defined multiplicatione routine
 * \param void (*Mult)(double* R,double* X,int n): pointer to system matrix multiplication method
 * 
 * There is a default matrix-multiplication routine implemented in the code but the user might want to 
 * use its own, wich can be set with this function. If the preconditioner mode is also "USER_DEFINED", then 
 * no explicit system matrix must be set. 
 * 
 * In the multiplicatione method the result is stored at pointer R and X represents the pointer to the vector that 
 * is multiplied by the system matrix. n denotes the dimension of the system.
 */
void PCG_set_mult(void (*Mult)(double* R,double* X,int n)){
	user_def_mult = Mult;
}




/** 
 * \brief Solves a linear equation with the PCG method. 
 * \param double* X: pointer to the vector where the solution of A.X=B is stored. 
 * \param double* B: pointer to the vector representing the right-hand side if the equation A.X=B
 * \param int n: dimension of the system
 * \param double tol: tolerance for the residuum reductioen r/r0 
 * \param int max_iter: maximum number of iterations until termination
 * \return int: actual number of iterations for the PCG to reduce the residuum ratio to r/r0<tol
 * 
 * This is the solution method for the linear system A.X=B with help of the Preconditioned Conjugated Gradient (PCG) 
 * method. Before using this routine, one must set the preconditioner and the system matrix. Both can be done 
 * by calling "PCG_set_preconditioner_mode". Alternatively, one can specify the matrix-mutiplcation and preconditioning 
 * routine by hand if favored. In the latter case, one does not have to pass the system matrix explicitly. 
 * 
 * When solving an equation, an initial guess X0 must be pointed to by X. After the algorithm has converged the 
 * solution is stored at X, and the initial guess is overwritten. Convergence is assumed when the ratio of the 
 * current residuum r/r0 = ||B-A.X||/||B-A.X0||<tol. If this condition is not fulfilled within max_iter iterations 
 * no convergence is assumed and the algorithm is terminated. 
 * 
 * Note that the PCG algorithm works only for symmetric and positive-definite matrices. The "ICHOL" 
 * preconditioner might complain for negative-definite matrices. In this case (not for indefinite matrices), 
 * simply multiply the matrix by -1 (i.e. solve (-A).x=-B). 
 */
int PCG_solve(double* X,double* B,int n,double tol,int max_iter){
	
	double r,alpha,beta,s,s_prev;
	double *P,*Q;
	
	double* R = zero_vector(n);
	(*PCG_mult)(R,X,n);
	vector_scale(R,-1.,n);
	vector_addmult(R,B,1.,n);
	
	double r0 = euklid_norm(R,n);
	if (r0==0){
		free(R);
		return 0;
	}
	
	double* Z = zero_vector(n);
	
	int i = 0;
	do{
		(*precon)(Z,R,n);
		s_prev = s;
		s = vector_scalar(Z,R,n);
		
		if (i>0){
			beta = s/s_prev;
			vector_scale(P,beta,n);
			vector_addmult(P,Z,1.,n);
		}
		else{
			P = copy_vector(Z,n);
			Q = zero_vector(n);
		}
		(*PCG_mult)(Q,P,n);
		alpha = s/vector_scalar(P,Q,n);
		vector_addmult(X,P,alpha,n);
		vector_addmult(R,Q,-alpha,n);		
		r = euklid_norm(R,n);
		i++;
		
		//printf("residuum ratio: %e\n",r/r0);	
		if (fabs(r/r0)<tol) break;
			
				
	}while(i<max_iter);
	pcg_residuum = r;
	
	free(P);
	free(Q);
	free(R);
	free(Z);	
	return i;
}
