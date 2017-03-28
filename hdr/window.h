#include <stdio.h>
#include <math.h>
#include <complex.h>

void window(unsigned int len, double(*win)(unsigned int, unsigned int), double complex * z) {
	if (len > 1) {
		int i;
		for (i = 0; i < len; i++) {
			z[i] *= win(i, len-1);
		}
	}
}

double blackman(unsigned int i, unsigned int len) {
	return 0.42 - 0.5*cos(2*i*M_PI/len) + 0.08*cos(4*i*M_PI/len);
}

double hamming(unsigned int i, unsigned int len) {
	return 0.54 - 0.46*cos(2*i*M_PI/len);
}

double hanning(unsigned int i, unsigned int len) {
	return 0.5 - 0.5*cos(2*i*M_PI/len);
}