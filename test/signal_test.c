#include <stdio.h>
#include <assert.h>
#include "../h/signal.h"

int main() {
	Segments *s = segments(0.9, /* PATH TO AUDIO FILE */, 4096);
	assert (NULL != s);
	printf("%d\n", s->segs[s->num-1][s->seglen-1]);
	
	return SUCCESS;
}