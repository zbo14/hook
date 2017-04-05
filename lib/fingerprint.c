#include <stdio.h>
#include <string.h>
#include <regex.h>
#include <sys/stat.h>
#include <unistd.h>
#include "../h/analysis.h"

int main(void) {
	printf("Enter a path to an audio file: ");
	char chr, *path;
	if (NULL == (path = malloc(sizeof(char)))) {
		return FAILURE;
	}
	path[0] = getchar();
	int i;
	for (i = 1; ;i++) {
		if ('\n' == (chr = getchar())) {
			break;
		}
		if (NULL == (path = realloc(path, (i+1)*sizeof(char)))) {
			return FAILURE;
		}
		path[i] = chr;
	}	
	Segments *s;
	if (NULL == (s = segments(OVERLAP_RATIO, path, WINDOW_SIZE))) {
		return FAILURE;
	}
	Egram *egram;
	if (NULL == (egram = default_entropygram(s))) {
		return FAILURE;
	}
	free(s);
	Fprint *fprint;
	if (NULL == (fprint = ses(egram))) {
		return FAILURE;
	}
	free(egram);
	regex_t regex;
	if (regcomp(&regex, "/([a-zA-Z0-9_-]*).mp3|wav$", REG_EXTENDED)) {
		return FAILURE;
	}
	regmatch_t match[2];
	if (regexec(&regex, path, 2, match, 0)) {
		return FAILURE;
	}
	regfree(&regex);
	char *filename;
	if (NULL == (filename = malloc((match[1].rm_eo-match[1].rm_so)*sizeof(char)))) {
		return FAILURE;
	}
	for (i = 0; i < match[1].rm_eo-match[1].rm_so; i++) {
		filename[i] = path[match[1].rm_so + i];
	}
	if (NULL == (path = realloc(path, 24 + strlen(filename)))) { // "output/" (7) + filename + "/fingerprint.bin" (16) + "\0" (1)
		return FAILURE;
	}
	strcpy(path, "output/");
	strcat(path, filename);
	mkdir(path, 0777);
	strcat(path, "/fingerprint.bin");
	puts(path);
	FILE *fp;
	if (NULL == (fp = fopen(path, "wb"))) {
		exit(FAILURE);
	}
	free(path);
	fwrite(&fprint->num, sizeof(unsigned int), 1, fp);
	for (i = 0; i < fprint->num; i++) {
		fwrite(fprint->segs[i], sizeof(uint8_t), SES_SIZE, fp);
	}
	fclose(fp);
	free(fprint);
	return SUCCESS;
}