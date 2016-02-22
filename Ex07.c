#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

void compute_v(int n, double *v) {

	int i;
	for ( i = 1; i < n+1; i++ ) {
		v[i-1] = 1/(double)( i*i );
	}
}

double compute_S(int n, double *v){

	double S = 0.0;
	int i;
	for ( i = 0; i < n; i++ ) {
		S += v[i];
	}
	return S;
}

double compute_error(double S_n) {

	double S = ( M_PI * M_PI )/6;
	return fabs(S - S_n);
}



int main(int argc, char ** argv) {


	MPI_Init(&argc, &argv);
	int P, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &P); // size = number of processes nprocs
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Numbering of process


	if (argc > 1) {
		// Only one process needs to print usage output
		if (rank == 0) {
			printf("Argument k needed.\n");
		}
		MPI_Finalize ();
		return 1;
	}

	double start = MPI_Wtime();
	sleep(2);

	int k = 14;
	int n = (int) pow(2, k);
	int np = (int) (n - n%P)/P;

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

		double err = compute_error(S);
		printf("Error = %1.16f\n", err);

		free(v);

	}
	else {

		double *v_n;
		if (rank == P-1) {
			v_n = calloc(np + n%P, sizeof(double));
		}
		else {
			v_n = calloc(np, sizeof(double));
		}
		MPI_Recv(v_n, np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

		double S_n = compute_S(np, v_n);

		free(v_n);

		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}
	
	double end = MPI_Wtime();
	
	if (rank == 0)
		printf("Elapsed time: %.2f seconds\n", end-start);
	
	
	MPI_Finalize();
	

	return 0;


}
