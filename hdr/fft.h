#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <complex.h>
#include "signal.h"


int fft_decimation_freq(unsigned int const len, double complex * z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int a, b, i, j, k;
	double complex cmplx1, cmplx2, temp1, temp2;
	for (a = len; a >= 2; a >>= 1) {
		b = a >> 1;
		cmplx1 = 1;
		cmplx2 = cos(M_PI/b) - sin(M_PI/b)*I;
		for (i = 0; i < b; i++) {
			for (j = i; j < len; j += a) {
				temp1 = z[j], temp2 = z[j+b];
				z[j] = temp1+temp2, z[j+b] = cmplx1*(temp1-temp2);
			}
			cmplx1 *= cmplx2;
		}
	}
	for (i = 1, j = 1; i <= len; i++) {
		if (i < j) {
			temp1 = z[i-1], temp2 = z[j-1];
			z[i-1] = temp2, z[j-1] = temp1;
		}
		for (k = len >> 1; k < j && k > 0; k >>= 1) {
			j -= k;
		}
		j += k;
	}
	return SUCCESS;
}


int inv_decimation_freq(unsigned int const len, double complex * z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int i;
	for (i = 0; i < len; i++) {
		z[i] = conj(z[i]);
	}
	if (FAILURE == fft_decimation_freq(len, z)) {
		return FAILURE;
	}
	for (i = 0; i < len; i++) {
		z[i] = conj(z[i]) / len;
	}
	return SUCCESS;
}

int fft_decimation_time(unsigned int const len, double complex *z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int a, b, i, j, k;
	int log2 = log2_floor(len);
	unsigned int sz = len * sizeof(double complex);
	double complex *q;
	if (NULL == (q = malloc(sz))) {
		return FAILURE;
	}
	memcpy(q, z, sz);
	for (i = 0; i < len; i++) {
		a = 0;
		for (j = 0; j < log2; j++) {
			a |= ((i >> (log2-j-1)) << j) & (1 << j);
		}
		z[i] = q[a%len];
	}
	free(q);
	double complex cmplx1, cmplx2, temp1, temp2;
	for (a = 2; a <= len; a <<= 1) {
		b = a >> 1;
		cmplx1 = 1;
		cmplx2 = cos(M_PI/b)-sin(M_PI/b)*I;
		for (i = 0; i < b; i++) {
			for (j = i; j < len; j += a) {
				temp1 = z[j], temp2 = cmplx1*z[j+b];
				z[j] = temp1+temp2, z[j+b] = temp1-temp2;
			}
		}
		cmplx1 *= cmplx2;
	}
	return SUCCESS;
}

int inv_decimation_time(unsigned int const len, double complex * z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int i;
	for (i = 0; i < len; i++) {
		z[i] = conj(z[i]);
	}
	if (FAILURE == fft_decimation_time(len, z)) {
		return FAILURE;
	}
	for (i = 0; i < len; i++) {
		z[i] = conj(z[i])/len;
	}
	return SUCCESS;
}

int fft_real(unsigned int const len, double complex *z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int len2 = len >> 1, len4 = len >> 2;
	unsigned int const sz = len * sizeof(double complex);
	double complex *p;
	if (NULL == (p = malloc(sz))) {
		return FAILURE;
	}
	int i;
	memcpy(p, z, sz);
	for (i = 0; i < len2; i++) {
		z[i] = p[i<<1] + p[(i<<1)+1]*I;
	}
	free(p);
	if (FAILURE == fft_decimation_freq(len2, z)) {
		return FAILURE;
	}
	double real1, real2, imag1, imag2;
	for (i = 1; i < len4; i++) {
		real1 = creal(z[i]), imag1 = cimag(z[i]);
		real2 = creal(z[len2-i]), imag2 = cimag(z[len2-i]);
		z[i+len2] = (imag1+imag2)/2 - (real1-real2)*I/2;
		z[len-i] = conj(z[i+len2]);
		z[i] = (real1+real2)/2 + (imag1-imag2)*I/2;
		z[len2-i] = conj(z[i]);
	}
	z[len*3/4] = cimag(z[len4]);
	z[len2] = cimag(z[0]);
	z[len4] = creal(z[len4]);
	z[0] = creal(z[0]);
	double complex cmplx1, cmplx2, temp1, temp2;
	cmplx1 = 1;
	cmplx2 = cos(M_PI/len2) - sin(M_PI/len2)*I;
	for (i = 0; i < len2; i++) {
		temp1 = z[i], temp2 = cmplx1*z[i+len2];
		z[i] = temp1+temp2, z[i+len2] = temp1-temp2;
		cmplx1 *= cmplx2;
	}
	return SUCCESS;
}

int inv_real(unsigned int const len, double complex *z) {
	if (!pow2(len)) {
		return FAILURE;
	}
	int i, len2 = len >> 1;
	for (i = 1; i < len2; i++) {
		z[i+len2] = conj(z[len2-i]);
	}
	for (i = 0; i < len; i++) {
		z[i] = creal(z[i]) + cimag(z[i]);
	}
	if (FAILURE == fft_real(len, z)) {
		return FAILURE;
	}
	for (i = 0; i < len; i++) {
		z[i] = (creal(z[i]) + cimag(z[i]))/len;
	}
	return SUCCESS;
}
