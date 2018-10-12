#include <iostream>

#include "qual_codec.h"


void fqz::encode_qual(RangeCoder *rc, char *qual, int len) {
    unsigned int last = 0;
    int delta = 5;
    int i, len2 = len;
    int q1 = 0, q2 = 0;

    /* Removing "Killer Bees" */
    while (len2 > 0 && qual[len2-1] == '#')
	len2--;

    for (i = 0; i < len2; i++) {
		unsigned char q = (qual[i] - '!') & (QMAX-1);

		model_qual[last].encodeSymbol(rc, q);

		// previous 2-3 bytes
		if (QBITS == 12) {
			last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
		} else {
			last = ((last << 6) + q) & ((1<<QBITS)-1);
		}

		if (qlevel > 1) {
			last  += (q1==q2) << QBITS;
			// delta saves 3-4%, but adds 14% cpu
			delta += (q1>q)*(q1-q);
			last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
		}

		if (qlevel > 2)
			//last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
			last += (MIN(i+15,127)&(15<<3))<<(QBITS+1);     // i>>3
			//last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

		_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
		q2 = q1; q1 = q;

		assert(last < (QSIZE*16));
    }

    if (len != len2)
	model_qual[last].encodeSymbol(rc, QMAX-1); /* terminator */
}

void fqz::decode_qual(RangeCoder *rc, char *qual, int len) {
    unsigned int last = 0;
    int i;
    int delta = 5;
    int q1 = 0, q2 = 0;

    for (i = 0; i < len; i++) {
		unsigned char q = model_qual[last].decodeSymbol(rc);

		if (q == QMAX-1) {
			while (i < len)
			qual[i++] = '#';
		} else {
			qual[i] = q + '!';

			if (QBITS == 12) {
				last = ((MAX(q1, q2)<<6) + q) & ((1<<QBITS)-1);
			} else {
				last = ((last << 6) + q) & ((1<<QBITS)-1);
			}

			if (qlevel > 1) {
				last  += (q1==q2) << QBITS;
				delta += (q1>q)*(q1-q);
				last  += (MIN(7*8, delta)&0xf8) << (QBITS-2);
			}

			if (qlevel > 2)
				//last += (MIN(i+7,127)&(7<<4))<<(QBITS);   // i>>4
				last += (MIN(i+15,127)&(15<<3))<<(QBITS+1);     // i>>3
				//last += (MIN(i+31,127)&(31<<2))<<(QBITS+2); // i>>2

			_mm_prefetch((const char *)&model_qual[last], _MM_HINT_T0);
			q2 = q1; q1 = q;
		}
	}
}
