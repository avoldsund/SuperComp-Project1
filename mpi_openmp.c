#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <stdbool.h>
#include <omp.h>

bool isPowerOfTwo(unsigned int x) {

    while( ((x % 2) == 0) && x > 1 ) {
        x /= 2;
    }

    return (x == 1);
}


int main(int argc, char ** argv) {

	MPI_Init(&argc, &argv);
	int P, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &P);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (isPowerOfTwo(P) == false) {
        if (rank == 0) {
            printf("The number of processes is not a power of two \n");
        }        
        MPI_Finalize();
        return 1;
    }

	if (argc < 2) {
		if (rank == 0) {
			printf("Usage: ex3 k, n=pow(2,k)\n");
			printf("n = vector/matrix size\n");
		}
		MPI_Finalize ();
		return 1;
	}

	double start = MPI_Wtime();

	
	int k = atoi(argv[1]);
	int n = (int) pow(2, k);
	
	int np = (int) (n - n%P)/P;
	
	double *v, S = 0.0;

	MPI_Status status;

	if ( rank == 0 ) {
		v = calloc(n, sizeof(double));
		
		int i;
		#pragma omp parallel for schedule(static)
		for ( i = 0; i < n; i++ ) {
			v[i] = 1 / (double)(i+1) / (double)(i+1);
		}
			
		for (i = 1; i < P; i++ ){
			MPI_Send( &v[np*i], np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		double S_n = 0.0;

		#pragma omp parallel for schedule(static) reduction(+:S_n)
		for ( i = 0; i < n; i++ ) {
			S_n += v[i];
		}

		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		double S = ( M_PI * M_PI )/6;
		double err = fabs(S - S_n);
		printf("Error = %.10f\n", err);

	}
	else {

		double *v_n;
		double S_n = 0.0;
		int i;
		if (rank == P-1) {
			v_n = calloc(np + n%P, sizeof(double));
			MPI_Recv(v_n, np+n%P, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			#pragma omp parallel for schedule(static) reduction(+:S_n)
			for ( i = 0; i < np + n%P; i++ ) {
				S_n += v_n[i];
			}
		}
		else {
			v_n = calloc(np, sizeof(double));
			MPI_Recv(v_n, np, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);

			#pragma omp parallel for schedule(static) reduction(+:S_n)
			for ( i = 0; i < np; i++ ) {
				S_n += v_n[i];
			}
		}

		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	double end = MPI_Wtime();
	
	if (rank == 0)
		printf("Elapsed time: %.10f seconds\n", end-start);

	MPI_Finalize();

	return 0;
}
