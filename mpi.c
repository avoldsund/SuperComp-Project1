#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#include <stdbool.h>

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

	
	int k = atoi(argv[1]);
	int n = (int) pow(2, k);
	
	int np = (int) (n - n%P)/P;
	
	double *v, S = 0.0;

	MPI_Status status;

	if ( rank == 0 ) {
		v = calloc(n, sizeof(double));
		
		int i;
		for ( i = 0; i < n; i++ ) {
			v[i] = 1 / double(i+1) / double(i+1);
		}
			
		for (i = 1; i < P; i++ ){
			MPI_Send( &v[np*i], np, MPI_DOUBLE, i, 0, MPI_COMM_WORLD);
		}
		double S_n = 0.0;
		for ( i = 0; i < n; i++ ) {
			S_n += v[i];
		}

		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
		
		double S = ( M_PI * M_PI )/6;
		double err = fabs(S - S_n);
		printf("Error = %f\n", err);

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
		MPI_Reduce(&S_n, &S, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

	return 0;
}