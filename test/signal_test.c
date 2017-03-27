#include <stdio.h>
#include <assert.h>
#include "../hdr/signal.h"

int main() {
	Segments *s = segments(0.9, "/Users/zach/Desktop/music/hey_jude_1.mp3", 4096);
	assert (NULL != s);
	printf("%d\n", s->segs[s->num-1][s->seglen-1]);
	
	return SUCCESS;
}