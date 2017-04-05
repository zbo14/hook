#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <complex.h>
#include "util.h"

FILE * get_file(char *path) {
	FILE *fp = fopen(path, "rb");
	if (NULL == fp) {
		return NULL;
	} 
	return fp;
}

unsigned int const file_length_16(FILE *fp) {
	fseek(fp, 0, SEEK_END);
	unsigned int const len = ftell(fp) / sizeof(uint16_t);
	rewind(fp);
	return len;
}

void read_time_domain_16(FILE *fp, unsigned int const len, unsigned int *x) {
	char hi, lo;
	int i;
	for (i = 0; i < len; i++) {
		hi = fgetc(fp), lo = fgetc(fp);
		if (feof(fp)) {
			break;
		}
		x[i] = uint16(hi, lo);
	}
	fclose(fp);
}

unsigned int const file_length_32(FILE *fp) {
	fseek(fp, 0, SEEK_END);
	unsigned int const len = ftell(fp) / sizeof(uint32_t);
	rewind(fp);
	return len;
}

void read_time_domain_32(FILE *fp, unsigned int const len, unsigned int *x) {
	char one, two, three, four;
	int i;
	for (i = 0; i < len; i++) {
		one = fgetc(fp), two = fgetc(fp), three = fgetc(fp), four = fgetc(fp);
		if (feof(fp)) {
			break;
		}
		x[i] = uint32(one, two, three, four);
	}
	fclose(fp);
}

unsigned int const num_segments(unsigned int const len, double olap, unsigned int const seglen) {
	if (olap < 0 || olap >= 1) {
		return FAILURE;
	}
	unsigned int numsegs;
	if (len == seglen) {
		numsegs = 1;
	} else if (len > seglen) {
		numsegs = (len-seglen)/(seglen-olap*seglen) + 2;
	} else {
		return FAILURE;
	}
	return numsegs;
}	

typedef struct _Segments32 {
	unsigned int len;
	unsigned int num;
	unsigned int seglen;
	unsigned int **segs;
} Segments;

Segments* segments(double olap, char *path, unsigned int const seglen, unsigned int size){
	if (size != 16 && size != 32) {
		return NULL;
	}
	FILE *fp;
	if (NULL == (fp = get_file(path))) {
		return NULL;
	}
	Segments *s;
	if (NULL == (s = malloc(sizeof(Segments)))) {
		return NULL;
	}
	unsigned int *x;
	if (NULL == (x = malloc(s->len * sizeof(unsigned int)))) {
		return NULL;
	}
	if (size == 16) {
		s->len = file_length_16(fp);
		read_time_domain_16(fp, s->len, x);
	} 
	if (size == 32) {
		s->len = file_length_32(fp);
		read_time_domain_32(fp, s->len, x);
	}
	s->num = num_segments(s->len, olap, seglen);
	s->seglen = seglen;
	if (NULL == (s->segs = malloc(s->num * sizeof(unsigned int*)))) {
		return NULL;
	}
	int i, j;
	for (i = 0, j = 0; i < s->num; i++, j += seglen-olap*seglen) {
		s->segs[i] = malloc(seglen * sizeof(unsigned int));
		memcpy(s->segs[i], x+j, seglen * sizeof(unsigned int));
	}
	free(x);
	return s;
}
