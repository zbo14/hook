#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#define SUCCESS 0
#define FAILURE -1

int pow2(int x) {
	return x != 0 && (x&(x-1))==0;
}

double square(double x) {
	return x*x;
}

int log2_floor(int x) {
	int i, log2 = 0;
	for (i = x; i > 1; i >>= 1) {
		log2++;
	} 
	return log2;
}

uint16_t uint16(char hi, char lo) {
	return (hi << 8) | lo;
}

uint32_t uint32(char one, char two, char three, char four) {
	return (one << 24) | (two << 16) | (three << 8) | four;
}

int* count_occurs(size_t len, uint16_t *x) {
	int *counts = malloc(len*sizeof(int));
	char *seen = malloc(len*sizeof(char));
	int i, j;
	for (i = 0; i < len; i++) {
		if (0 == seen[i]) {
            counts[i] = 1;
			for (j = i+1; j < len; j++) {
				if (x[i] == x[j]) {
					seen[j] = 1;
					counts[i]++;
				}
			}
		}
	}
	free(seen);
    return counts;
}
