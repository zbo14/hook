#include <stdio.h>
#include <sys/stat.h>
#include "../hdr/util.h"

int main(void) {
	if (0 != mkdir("output", 0777)) {
		perror("mkdir");
		return FAILURE;
	}
	return SUCCESS;
}