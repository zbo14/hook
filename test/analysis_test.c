#include <stdio.h>
#include <assert.h>
#include "../h/analysis.h"

int main() {
	size_t const len = 3;
	uint8_t b1[len] = {1, 2, 3}, b2[len] = {4, 5, 6};
	assert (7 == hamming_dist(b1, b2, len));
	
	Segments *s = segments(0.9, /* PATH TO AUDIO FILE */, 8192);
	Egram *egram = default_entropygram(s);
	Fprint *fprint1 = ses(egram);

	s = segments(0.9, /* PATH TO ANOTHER AUDIO FILE */, 8192);
	egram = default_entropygram(s);
	Fprint *fprint2 = ses(egram);

	double diff = dtw(fprint1, fprint2);
	printf("DTW DIFF: %f\n", diff);

	double sim = default_lcs(fprint1, fprint2);
	printf("LCS SIM: %f\n", sim);

	diff = default_levenshtein(fprint1, fprint2);
	printf("LEVENSHTEIN DIFF: %f\n", diff);

	free(fprint1);
	free(fprint2);

	return SUCCESS;
}