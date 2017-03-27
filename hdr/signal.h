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

size_t const file_length(FILE *fp) {
	fseek(fp, 0, SEEK_END);
	size_t const len = ftell(fp) / sizeof(uint16_t);
	rewind(fp);
	return len;
}

void read_time_domain(FILE *fp, size_t const len, uint16_t *x) {
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

size_t const num_segments(size_t const len, double olap, size_t const seglen) {
	if (olap < 0 || olap >= 1) {
		return FAILURE;
	}
	size_t numsegs;
	if (len == seglen) {
		numsegs = 1;
	} else if (len > seglen) {
		numsegs = (len-seglen)/(seglen-olap*seglen) + 2;
	} else {
		return FAILURE;
	}
	return numsegs;
}	

typedef struct _Segments {
	size_t len;
	size_t num;
	size_t seglen;
	uint16_t **segs;
} Segments;

Segments* segments(double olap, char *path, size_t const seglen){
	FILE *fp;
	if (NULL == (fp = get_file(path))) {
		return NULL;
	}
	Segments *s;
	if (NULL == (s = malloc(sizeof(Segments)))) {
		return NULL;
	}
	s->len = file_length(fp);
	uint16_t *x;
	if (NULL == (x = malloc(s->len * sizeof(uint16_t)))) {
		return NULL;
	}
	read_time_domain(fp, s->len, x);
	s->num = num_segments(s->len, olap, seglen);
	s->seglen = seglen;
	if (NULL == (s->segs = malloc(s->num * sizeof(uint16_t*)))) {
		return NULL;
	}
	int i, j;
	for (i = 0, j = 0; i < s->num; i++, j += seglen-olap*seglen) {
		s->segs[i] = malloc(seglen * sizeof(uint16_t));
		memcpy(s->segs[i], x+j, seglen * sizeof(uint16_t));
	}
	free(x);
	return s;
}
