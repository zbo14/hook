#include <stdio.h>
#include <assert.h>
#include <complex.h>
#include "../h/fft.h"

int main() {
	size_t const len = 16;
	double complex *z;
	if (NULL == (z = malloc(len * sizeof(double complex)))) {
		exit(FAILURE);
	}
	int i;
	for (i = 0; i < len; i++) {
		z[i] = (double)i + (double)i/len;
	}
	assert(SUCCESS == fft_real(len, z));
	assert(SUCCESS == inv_real(len, z));
	for (i = 0; i < len; i++) {
		printf("%f == %f\n", (double)i + (double)i/len, creal(z[i]));
	}

	return SUCCESS;
}