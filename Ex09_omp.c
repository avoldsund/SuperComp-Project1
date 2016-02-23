#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
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


int main(int argc, char **argv) {
	
	double start = omp_get_wtime();
	
	int n, k = atoi(argv[1]);
	
		
	n = (int) pow(2, k);

	double *v;
	v = calloc(n, sizeof(double));
	compute_v(n, v); 

	double S_n = compute_S(n, v);

	free(v);

	double err = compute_error(S_n);
	printf("Error = %.16f\n", err);

	double end = omp_get_wtime();
	printf("Elapsed time: %.8f seconds\n", end-start);


	return 0;


}
