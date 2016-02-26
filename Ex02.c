#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void generate_v(int n, double *v) {

	int i;	
	for ( i = 0; i < n; i++ ) {
		v[i] = 1/ (i+1) / (i+1) ;
	}
}

double compute_S(int n, double *v){

	double S = 0.0;
	int i;

	//#pragma omp parallel for schedule(static) reduction(+:S)
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
		printf("err[k = %d] = %.16f\n", i+3, vec[i]);
	}
}


int main() {

	double start, end;
	start = omp_get_wtime();
	
	int n, k;
	//int k_max = 15, k_min = 3;

	//double err_vec[k_max - k_min];
	
	//for ( k = k_min; k < k_max; k++ ) {
		
 		n = (int) pow(2, k);
		
		double *v;
		v = calloc(n, sizeof(double));
		generate_v(n, v); 

		double S_n = compute_S(n, v);

		free(v);

		double err = compute_error(S_n);
		//err_vec[k-3] = err;
	//}
	end = omp_get_wtime();

	//print_vec(k_max - k_min, err_vec);
	printf("Error = %.10f\n", err);
	printf("Elapsed time: %.10f\n", end-start);

	return 0;

}
