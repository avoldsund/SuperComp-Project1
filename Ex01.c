#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void generate_v(int n, double *v) {

	int i;	

	for ( i = 0; i < n; i++ ) {
		v[i] = 1/(double)( (i+1)*(i+1) );
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

void print_vec(int n, double *vec) {
	
	int i;

	for ( i = 0; i < n; i++ ) {
		printf("err[k = %d] = \t %1.16f\n", i+3, vec[i]);
	}
}


int main() {
	
	int n, k;
	int k_max = 15, k_min = 3;

	double err_vec[k_max - k_min];
	
	for ( k = k_min; k < k_max; k++ ) {
		
 		n = (int) pow(2, k);
		
		double *v;
		v = calloc(n, sizeof(double));
		generate_v(n, v); 

		double S_n = compute_S(n, v);

		free(v);

		double err = compute_error(S_n);
		err_vec[k-3] = err;
	}

	print_vec(k_max - k_min, err_vec);


	return 0;
}
