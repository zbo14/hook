#include <stdio.h>
#include <assert.h>
#include "../hdr/analysis.h"

int main() {
	size_t const len = 3;
	uint8_t b1[len] = {1, 2, 3}, b2[len] = {4, 5, 6};
	assert (7 == hamming_dist(b1, b2, len));

	/*
	Segments *s = segments(0.9, "/Users/zach/Desktop/music/rhapsody_1.mp3", 4096);
	Sgram *sgram = default_spectrogram(s);
	Peaks *peaks = default_find_peaks(sgram);
	Constellation *c1 = default_constellate(peaks);

	s = segments(0.5, "/Users/zach/Desktop/music/rhapsody_2.mp3", 4096);
	sgram = default_spectrogram(s);
	peaks = default_find_peaks(sgram);
	Constellation *c2 = default_constellate(peaks);

	double score = compare_constellations(c1, c2, 0.1);
	printf("score: %f\n", score);
	*/

	Segments *s = segments(0.9, "/Users/zach/Desktop/music/rhapsody_1.mp3", 8192);
	Egram *egram = default_entropygram(s);
	EFprint *fprint1 = efp(egram);

	s = segments(0.9, "/Users/zach/Desktop/music/rhapsody_2.mp3", 8192);
	egram = default_entropygram(s);
	EFprint *fprint2 = efp(egram);
	
	double diff = dtw(fprint1, fprint2);
	printf("DTW DIFF: %f\n", diff);

	double sim = lcs(0.2, fprint1, fprint2);
	printf("LCS SIM: %f\n", sim);

	diff = levenshtein(0.2, fprint1, fprint2);
	printf("LEVENSHTEIN DIFF: %f\n", diff);

	return SUCCESS;
}