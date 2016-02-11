#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <omp.h>

void compute_v(int n, double *v) {
	int i;
	
	for ( i = 0; i < n; i++ ) {
		v[i] = 1/(double)( (i+1)*(i+1) );
	}
}

double compute_S(int n, double *v){

	double S = 0.0;
	int i;

	#pragma omp parallel for schedule(static) reduction(+:S)
	for ( i = 0; i < n; i++ ) {
		S += v[i];
	}
	return S;
}

double compute_error(double S_n) {

	double S = ( M_PI * M_PI )/6;
	return fabs(S - S_n);
}

void print_vec(int n, double *vec) {
	
	int i;
	for ( i = 0; i < n; i++ ) {
		printf("err[k = %d] = %f\n", i+3, vec[i]);
	}
}


int main(int argc, char ** argv) {


	MPI_Init(&argc, &argv);
	int P, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &P); // size = number of processes nprocs
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Numbering of process

	if (argc < 2) {
		// Only one process needs to print usage output
		if (rank == 0) {
			printf("Usage: ex3 k, n=pow(2,k)\n");
			printf("n = vector/matrix size\n");
		}
		MPI_Finalize ();
		return 1;
	}

	
	int k = atoi(argv[1]);
	int n = (int) pow(2, k);
	
	int np = (int) (n - n%P)/P;
	//printf("n mod P = %d\n", n%P);
	
	double *v, S = 0.0;

	MPI_Status status;
	int i;	

	if ( rank == 0 ) {
		v = calloc(n, sizeof(double));
		compute_v(n, v);
			
		for (i = 1; i < P; i++ ){
			MPI_Send( &v[np*i], np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		double S_n = compute_S(np, v);
		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		printf("S = %f\n", S);
		double err = compute_error(S);
		printf("Error = %f\n", err);

		//int n_check = (P)*np + n%P;

	}
	else {

		double *v_n;
		if (rank == P-1) {
			v_n = calloc(np + n%P, sizeof(double));
			//printf("np + mod = %d\n", np + n%P);
		}
		else {
			v_n = calloc(np, sizeof(double));
			//printf("np = %d\n", np);
		}
		MPI_Recv(v_n, np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);



		double S_n = compute_S(np, v_n);
		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	
	MPI_Finalize();


	return 0;


}
