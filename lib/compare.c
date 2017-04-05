#include <stdio.h>
#include "../hdr/analysis.h"

int main(void) {
	printf("Enter name of fingerprinted recording: ");
	char chr, *name1, *path1;
	if (NULL == (name1 = malloc(sizeof(char)))) {
		return FAILURE;
	}
	name1[0] = getchar();
	int i, j;
	for (i = 1; ;i++) {
		if ('\n' == (chr = getchar())) {
			break;
		}
		if (NULL == (name1 = realloc(name1, (i+1)*sizeof(char)))) {
			return FAILURE;
		}
		name1[i] = chr;
	}
	if (NULL == (path1 = malloc(24 + strlen(name1)))) {
		return FAILURE;
	}
	strcpy(path1, "output/");
	strcat(path1, name1);
	strcat(path1, "/fingerprint.bin");
	char *name2, *path2;
	printf("Enter name of another fingerprinted recording: ");	
	if (NULL == (name2 = malloc(sizeof(char)))) {
		return FAILURE;
	}
	name2[0] = getchar();
	for (i = 1; ;i++) {
		if ('\n' == (chr = getchar())) {
			break;
		}
		if (NULL == (name2 = realloc(name2, (i+1)*sizeof(char)))) {
			return FAILURE;
		}
		name2[i] = chr;
	}
	if (NULL == (path2 = malloc(24 + strlen(name2)))) {
		return FAILURE;
	}
	strcpy(path2, "output/");
	strcat(path2, name2);
	strcat(path2, "/fingerprint.bin");
	FILE *fp1, *fp2;
	if (NULL == (fp1 = fopen(path1, "rb"))) {
		return FAILURE;
	}
	if (NULL == (fp2 = fopen(path2, "rb"))) {
		return FAILURE;
	}
	Fprint *fprint1, *fprint2;
	if (NULL == (fprint1 = malloc(sizeof(Fprint)))) {
		return FAILURE;
	}
	fread(&fprint1->num, sizeof(unsigned int), 1, fp1);
	if (NULL == (fprint1->segs = malloc(fprint1->num * sizeof(uint8_t*)))) {
		return FAILURE;
	}
	for (i = 0; i < fprint1->num; i++) {
		if (NULL == (fprint1->segs[i] = malloc(SES_SIZE * sizeof(uint8_t)))) {
			return FAILURE;
		}
		fread(fprint1->segs[i], sizeof(uint8_t), SES_SIZE, fp1);
	}
	for (i = 0; i < fprint1->num; i++) {
		for (j = 0; j < SES_SIZE; j++) {
			printf("%d\n", fprint1->segs[i][j]);
		}
		printf("\n");
	}
	fclose(fp1);
	if (NULL == (fprint2 = malloc(sizeof(Fprint)))) {
		return FAILURE;
	}
	fread(&fprint2->num, sizeof(unsigned int), 1, fp2);
	if (NULL == (fprint2->segs = malloc(fprint2->num * sizeof(uint8_t*)))) {
		return FAILURE;
	}
	for (i = 0; i < fprint2->num; i++) {
		if (NULL == (fprint2->segs[i] = malloc(SES_SIZE * sizeof(uint8_t)))) {
			return FAILURE;
		}
		fread(fprint1->segs[i], sizeof(uint8_t), SES_SIZE, fp1);
	}
	fclose(fp2);

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