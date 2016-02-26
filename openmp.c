#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>


int main(int argc, char ** argv) {

	if (argc < 2) {
		printf("Usage: k, n=pow(2,k)\n");
		printf("n = vector/matrix size\n");
		return 1;
	}


	double start;
	start = omp_get_wtime();
	sleep(2);
	
	int k = atoi(argv[1]);	
	int n = (int) pow(2, k);
	
	double *v;
	v = (double*) malloc(n * sizeof(double));

	int i;
	#pragma omp parallel for schedule(static)
	for ( i = 0; i < n; i++ ) {
		v[i] = 1/ (double)(i+1)/(double)(i+1);
	}

	double S_n = 0.0;
	//#pragma omp parallel for schedule(static) reduction(+:S_n)
	for ( i = 0; i < n; i++ ) {
		S_n += v[i];
	}

	free(v);

	double S = ( M_PI * M_PI )/6;
	double err = fabs(S - S_n);

	double end;
	end = omp_get_wtime();

	printf("Error = %.10f\n", err);
	printf("Elapsed time: %.10f\n", end-start-2.0);

	return 0;

}
