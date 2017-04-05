#include <stdio.h>
#include <assert.h>
#include "../h/util.h"

int main() {
	assert(3 == log2_floor(9));
	assert(6 == log2_floor(64));
	assert (pow2(8));
	assert (!pow2(65));
	assert (square(8) == 64);
	assert (square(10) == 100);
	assert (uint16(1, 1) == 257);
	assert (uint16(255, 255) == 65535);

	size_t const len = 10;
	uint16_t x[len] = {1,1,2,3,4,3,2,3,4,2};
    int *counts = count_occurs(len, x);
   	assert (2 == counts[0]);
   	assert (3 == counts[2]);
   	assert (3 == counts[3]);
   	assert (2 == counts[4]);

   	return SUCCESS;
}