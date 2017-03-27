#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "fft.h"
#include "window.h"

#define BANDS 24
#define EFP_SIZE 3
#define FAN_VALUE 15
#define NEIGHBORHOOD_SIZE 20
#define OVERLAP_RATIO 0.9
#define SAMPLING_RATE 44100
#define SFP_SIZE 24
#define TIME_DELTA 200
#define WINDOW_SIZE 8192

int hamming_dist(uint8_t *bytes1, uint8_t *bytes2, size_t len) {
	int dist = 0, i, j;
	for (i = 0; i < len; i++) {
		for (j = 0; j < 8; j++) {
			if ((bytes1[i]&(1<<j))>>j != (bytes2[i]&(1<<j))>>j) {
				dist++;
			}
		}
	}
	return dist;
}

double freq_to_bark(double freq) {
	double bark = 26.81*freq/(1960+freq) - 0.53;
	if (bark < 2) {
		bark += 0.15 * (2 - bark);
	} else if (bark > 20.1) {
		bark += 0.22 * (bark - 20.1);
	}
	return bark;
}

typedef struct _Egram {
	size_t num;
	size_t seglen;
	double **segs;
} Egram;

Egram* entropygram(unsigned int fs, Segments *s, double(*win)(unsigned int, size_t)) {
	int seglen2 = (s->seglen)>>1;
	int *bands;
	if (NULL == (bands = malloc(BANDS * sizeof(int)))) {
		return NULL;
	}
	int band, i;
	for (band = 1, i = 0; band <= BANDS && i <= seglen2; i++) {
		if (band < freq_to_bark((double)i*fs/(s->seglen))) {
			bands[band-1] = i;
			band++;
		}
	}
	Egram *egram;
	if (NULL == (egram = malloc(sizeof(Egram)))) {
		return NULL;
	}
	egram->num = s->num;
	egram->seglen = s->seglen;
	if (NULL == (egram->segs = malloc(egram->num * sizeof(double*)))) {
		return NULL;
	}
	int end, j, k, start;
	double avgImag, avgReal, covar, varImag, varReal;
	double complex *z;
	if (NULL == (z = malloc(egram->seglen * sizeof(double complex)))) {
		return NULL;
	}
	for (i = 0; i < egram->num; i++) {
		for (j = 0; j < egram->seglen; j++) {
			z[j] = s->segs[i][j];
		}
		free(s->segs[i]);
		window(egram->seglen, win, z);
		if (FAILURE == fft_real(egram->seglen, z)) {
			return NULL;
		}
		if (NULL == (egram->segs[i] = malloc(BANDS * sizeof(double)))) {
			return NULL;
		}
		avgImag = 0, avgReal = 0, start = 0;
		for (j = 0; j < BANDS; j++) {
			end = bands[j];
			for (k = start; k < end; k++) {
				avgImag += cimag(z[k]) / (seglen2+1);
				avgReal += creal(z[k]) / (seglen2+1);
			}
			covar = 0, varImag = 0, varReal = 0;
			for (k = start; k < end; k++) {
				covar += (cimag(z[k])-avgImag) * (creal(z[k])-avgReal) / seglen2;
				varImag += square(cimag(z[k])-avgImag) / seglen2;
				varReal += square(creal(z[k])-avgReal) / seglen2;
			}
			egram->segs[i][j] = log(2*M_PI*M_E) + 0.5*log(varImag*varReal-square(covar));
			start = end;
		}
	}
	free(bands);
	free(s);
	free(z);
	return egram;
}


Egram* default_entropygram(Segments *s) {
	return entropygram(SAMPLING_RATE, s, hanning);
}

typedef struct _EFprint {
	size_t num;
	uint8_t **segs;
} EFprint;

EFprint* efp(Egram *egram) {
	EFprint *fprint;
	if (NULL == (fprint = malloc(sizeof(EFprint)))) {
		return NULL;
	}
	if (NULL == (fprint->segs = malloc(egram->num * sizeof(uint8_t*)))) {
		return NULL;
	}
	fprint->num = egram->num;
	int i, j;
	for (i = 0; i < fprint->num; i++) {
		fprint->segs[i] = malloc(EFP_SIZE * sizeof(uint8_t));
		if (i == 0) {
			for (j = 0; j < EFP_SIZE; j++) {
				fprint->segs[i][j] = 0;
			}
			continue;
		} 
		for (j = 0; j < BANDS; j++) {
			if (egram->segs[i][j] > egram->segs[i-1][j]) {
				fprint->segs[i][j/8] |= (1 << (j%8));
			}
		}
		if (i > 0) {
		 	free(egram->segs[i-1]);
		}
	}
	free(egram);
	return fprint;
}

double dtw(EFprint *fprint1, EFprint *fprint2) {
	int *col1, *col2;  
	if (NULL == (col1 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	if (NULL == (col2 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	int i, j;
	for (i = 0; i < fprint1->num; i++) {
		for (j = 0; j < fprint2->num; j++) {
			if (i == 0 || j == 0) {
				col2[j] = hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE);
			} else {
				if (col1[j-1] < col1[j] && col1[j-1] < col2[j-1]) {
					col2[j] = col1[j-1] + 2*hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE);
				} else if (col1[j] < col1[j-1] && col1[j] < col2[j-1]) {
					col2[j] = col1[j] + hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE);
				} else {
					col2[j] = col2[j-1] + hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE);
				}
			}
		}
		col1 = col2;
	}
	double diff = (double)col1[fprint2->num-1]/(fprint1->num+fprint2->num);
	free(col1);
	return diff;
}

double lcs(double eps, EFprint *fprint1, EFprint *fprint2) {
	int *col1, *col2;  
	if (NULL == (col1 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	if (NULL == (col2 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	int i, j;
	for (i = 0; i <= fprint2->num; i++) {
		col1[i] = i;
	}
	for (i = 0; i < fprint1->num; i++) {
		for (j = 0; j < fprint2->num; j++) {
			if (eps > hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE)) {
				col2[j+1] = col1[j] + 1;
			} else if (col2[j] > col1[j+1]) {
				col2[j+1] = col2[j];
			} else {
				col2[j+1] = col1[j+1];
			}
		}
		col1 = col2;
	}
	double sim = (double)col1[fprint2->num]/(fprint1->num+fprint2->num);
	free(col1);
	return sim;
}

double levenshtein(double eps, EFprint *fprint1, EFprint *fprint2) {
	int *col1, *col2;  
	if (NULL == (col1 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	if (NULL == (col2 = malloc(fprint2->num * sizeof(int)))) {
		return FAILURE;
	}
	int dist, i, j;
	for (i = 0; i <= fprint2->num; i++) {
		col1[i] = i;
	}
	for (i = 0; i < fprint1->num; i++) {
		for (j = 0; j < fprint2->num; j++) {
			if ((dist = hamming_dist(fprint1->segs[i], fprint2->segs[j], EFP_SIZE)), dist < eps) {
				col2[j+1] = col1[j];
			} else if (col1[j] < col2[j] && col1[j] < col1[j+1]) {
				col2[j+1] = col1[j] + dist;
			} else if (col2[j] < col1[j] && col2[j] < col1[j+1]) {
				col2[j+1] = col2[j] + 1;
			} else {
				col2[j+1] = col1[j+1] + 1;
			}
		}
		col1 = col2;
	}
	double diff;
	if (fprint1->num > fprint2->num) {
		diff = (double)col1[fprint2->num]/fprint1->num;
	} else {
		diff = (double)col1[fprint2->num]/fprint2->num;
	}
	free(col1);
	return diff;
}

//---------------------------------------------------------------

typedef struct _Sgram {
	double *freqs;
	size_t num;
	size_t seglen;
	double **segs;
} Sgram;

Sgram* spectrogram(unsigned int fs, Segments *s, double(*win)(unsigned int, size_t)) {
	int seglen2 = (s->seglen)>>1;
	Sgram* sgram;
	if (NULL == (sgram = malloc(sizeof(Sgram)))) {
		return NULL;
	}
	if (NULL == (sgram->freqs = malloc((seglen2+1)*sizeof(double)))) {
		return NULL;
	}	
	sgram->num = s->num;
	sgram->seglen = seglen2 + 1;
	if (NULL == (sgram->segs = malloc(sgram->num * sizeof(double*)))) {
		return NULL;
	}
	double complex *z;
	if (NULL == (z = malloc(s->seglen * sizeof(double complex)))) {
		return NULL;
	}
	int i, j; 
	for (i = 0; i < sgram->num; i++) {
		for (j = 0; j < s->seglen; j++) {
			z[j] = s->segs[i][j];
		}
		free(s->segs[i]);
		window(s->seglen, win, z);
		if (FAILURE == fft_real(s->seglen, z)) {
			return NULL;
		}
		if (NULL == (sgram->segs[i] = malloc(sgram->seglen*sizeof(double)))) {
			return NULL;
		}
		for (j = 0; j < sgram->seglen; j++) {
			sgram->freqs[j] = (double)j*fs/s->seglen;
			sgram->segs[i][j] = log(creal(z[j] * conj(z[j])));
		}
	}
	free(s);
	free(z);
	return sgram;
}

Sgram* default_spectrogram(Segments *s) {
	return spectrogram(SAMPLING_RATE, s, hanning);
}

typedef struct _Peaks {
	size_t num;
	size_t *seglens;
	double **segs;
} Peaks;

Peaks* find_peaks(unsigned int nbr, Sgram *sgram) {
	Peaks *peaks;
	if (NULL == (peaks = malloc(sizeof(Peaks)))) {
		return NULL;
	}
	peaks->num = sgram->num;
	if (NULL == (peaks->seglens = malloc(peaks->num * sizeof(size_t)))) {
		return NULL;
	}
	if (NULL == (peaks->segs = malloc(peaks->num * sizeof(double*)))) {
		return NULL;
	}
	int i, j, k, l, seglen;
	for (i = 0; i < sgram->num; i++) {
		seglen = 0;
		for (j = 0; j < sgram->seglen; ) {
			for (k = 0; k <= nbr; k++) {
				if (k > 0) {
					if (j-k >= 0) {
						if (sgram->segs[i][j] < sgram->segs[i][j-k]) {
							goto end;
						}
					}
					if (j+k < sgram->seglen) {
						if (sgram->segs[i][j] < sgram->segs[i][j+k]) {
							goto end;
						}
					}
				}
				for (l = 1; l <= nbr-k; l++) {
					if (i-l >= 0) {
						if (j-k >= 0) {
							if (sgram->segs[i][j] <= sgram->segs[i-l][j-k]) {
								goto end;
							}
						}
						if (j+k < sgram->seglen) {
							if (sgram->segs[i][j] <= sgram->segs[i-l][j+k]) {
								goto end;
							}
						}
					}
					if (i+l < sgram->num) {
						if (j-k >= 0) {
							if (sgram->segs[i][j] <= sgram->segs[i+l][j-k]) {
								goto end;
							}
						}
						if (j-k >= 0) {
							if (sgram->segs[i][j] <= sgram->segs[i+l][j+k]) {
								goto end;
							}
						}
					}
				}
			}
			if (seglen++, seglen == 1) {
				if (NULL == (peaks->segs[i] = malloc(sizeof(double)))) {
					return NULL;
				}
			} else {
				if (NULL == (peaks->segs[i] = realloc(peaks->segs[i], seglen*sizeof(double)))) {
					return NULL;
				}
			}
			peaks->segs[i][seglen-1] = sgram->freqs[j];
			end: j++;
		}
		peaks->seglens[i] = seglen;
	}
	free(sgram);
	return peaks;
}

Peaks* default_find_peaks(Sgram *sgram) {
	return find_peaks(NEIGHBORHOOD_SIZE, sgram);
}

// TODO: calc sha3 checksum

uint8_t* sfp(double peak1, double peak2, uint64_t tdelta) {
	uint8_t *fprint;
	if (NULL == (fprint = malloc(2*sizeof(double) + sizeof(uint64_t)))) {
		return NULL;
	}
	memcpy(fprint, &peak1, sizeof(double));
	memcpy(fprint + sizeof(double), &peak2, sizeof(double));
	memcpy(fprint + 2*sizeof(double), &tdelta, sizeof(uint64_t));
	return fprint;
}

typedef struct _Constellation {
	uint8_t **fprints;
	size_t num;
} Constellation;

Constellation* constellate(unsigned int fan, Peaks *peaks, uint64_t tdelta) {
	Constellation *c;
	if (NULL == (c = malloc(sizeof(Constellation)))) {
		return NULL;
	}
	c->num = 0;
	int count, i, j, k, l;
	for (i = 0; i < peaks->num; i++) {
		for (j = 0; j < peaks->seglens[i]; j++) {
			for (count = 1, k = j+1; k < peaks->seglens[i]; count++, k++) {
				if (c->num++, c->num == 1) {
					if (NULL == (c->fprints = malloc(sizeof(uint8_t*)))) {
						return NULL;
					}
				} else {
					if (NULL == (c->fprints = realloc(c->fprints, c->num*sizeof(uint8_t*)))) {
						return NULL;
					}
				}
				c->fprints[c->num-1] = sfp(peaks->segs[i][j], peaks->segs[i][k], 0);
				if (count == fan) {
					goto end;
				}
			}
			for (k = i+1; k < peaks->num && k-i <= tdelta; count++, k++) {
				for (l = 0; l < peaks->seglens[k]; l++) {
					c->num++;
					if (NULL == (c->fprints = realloc(c->fprints, c->num*sizeof(uint8_t*)))) {
						return NULL;
					}
					c->fprints[c->num-1] = sfp(peaks->segs[i][j], peaks->segs[k][l], k-i);
					if (count == fan) {
						goto end;
					}
				}
			}
			end:;
		}
	}
	return c;
}

Constellation* default_constellate(Peaks *peaks) {
	return constellate(FAN_VALUE, peaks, TIME_DELTA);
}

int twlcs(Constellation *c1, Constellation *c2) {
	int *col1, *col2; 
	if (NULL == (col1 = malloc((c2->num+1)*sizeof(int)))) {
		return FAILURE;
	}
	if (NULL == (col2 = malloc((c2->num+1)*sizeof(int)))) {
		return FAILURE;
	}
	int i, j, k;
	for (i = 0; i <= c2->num; i++) {
		col1[i] = i;
	}
	for (i = 0; i < c1->num; i++) {
		for (j = 0; j < c2->num; j++) {
			for (k = 0; k < SFP_SIZE; k++) {
				if (c1->fprints[i][k] != c2->fprints[j][k]) {
					if (col2[j] > col1[j+1]) {
						col2[j+1] = col2[j];
					} else {
						col2[j+1] = col1[j+1];
					}
					goto end;
				}
			}
			if (col1[j] > col2[j] && col1[j] > col1[j+1]) {
				col2[j+1] = col1[j] + 1;
			} else if (col2[j] > col1[j] && col2[j] > col1[j+1]) {
				col2[j+1] = col2[j] + 1;
			} else {
				col2[j+1] = col1[j+1] + 1;
			}
			end: col1 = col2;
		}
	}
	int seqlen = col1[c2->num];
	free(col1);
	return seqlen;
}

double compare_constellations(Constellation *c1, Constellation *c2, double eps) {
	double score;
	if (c1->num < c2->num) {
		score = c1->num;
	} else {
		score = c2->num;
	}
	int seqlen1, seqlen2;
	if (FAILURE == (seqlen1 = twlcs(c1, c2))) {
		return FAILURE;
	}
	if (FAILURE == (seqlen2 = twlcs(c2, c1))) {
		return FAILURE;
	}
	if (seqlen1 > seqlen2) {
		score = seqlen1 / score;
	} else {
		score = seqlen2 / score;
	}
	return score;
}
