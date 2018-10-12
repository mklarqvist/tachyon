#include <limits.h>
#include <cstdint>
#include <cassert>
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#define  DO(n) for (int _=0; _<n; _++)
#define  TOP   (1<<24)

#define ABS(a)   ((a)>0?(a):-(a))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

class RangeCoder {
private:
    uint64_t low;
    uint32_t range, code;

public:
	uint8_t *in_buf;
	uint8_t *out_buf;

	void input (char *in)  { out_buf = in_buf = (uint8_t*)in;  }
	void output(char *out) { in_buf = out_buf = (uint8_t*)out; }
	char* input(void)  {return (char *)in_buf;}
	char* output(void) {return (char *)out_buf;}
	int size_out(void) {return out_buf - in_buf;}
	int size_in(void)  {return in_buf - out_buf;}

	void StartEncode ( void )
	{
		low=0;
		range=(uint32_t)-1;
	}

	void StartDecode( void )
	{
		low=0;
		range=(uint32_t)-1;
		DO(8) code = (code<<8) | *in_buf++;
	}

	void FinishEncode( void ) { DO(8) (*out_buf++ = low>>56), low<<=8; }
	void FinishDecode( void ) {}

	void Encode (uint32_t cumFreq, uint32_t freq, uint32_t totFreq)
	{
		low  += cumFreq * (range/= totFreq);
		range*= freq;

		if (cumFreq + freq > totFreq)
			abort();

		while( range<TOP ) {
			// range = 0x00ffffff..
			// low/high may be matching
			//       eg 88332211/88342211 (range 00010000)
			// or differing
			//       eg 88ff2211/89002211 (range 00010000)
			//
			// If the latter, we need to reduce range down
			// such that high=88ffffff.
			// Eg. top-1      == 00ffffff
			//     low|top-1  == 88ffffff
			//     ...-low    == 0000ddee
			if ( uint8_t((low^(low+range))>>56) )
			range = ((uint32_t(low)|(TOP-1))-uint32_t(low));
			*out_buf++ = low>>56, range<<=8, low<<=8;
		}
	}

	uint32_t GetFreq (uint32_t totFreq) {
		return code/(range/=totFreq);
	}

	void Decode (uint32_t cumFreq, uint32_t freq, uint32_t totFreq)
	{
		uint32_t temp = cumFreq*range;
		low  += temp;
		code -= temp;
		range*= freq;

		while( range<TOP ) {
			if ( uint8_t((low^(low+range))>>56) )
				range = ((uint32_t(low)|(TOP-1))-uint32_t(low));
			code = (code<<8) | *in_buf++, range<<=8, low<<=8;
		}
	}
};

// Start simple model
// Shrinking this to 1<<10 gives 2-3% smaller qualities, but 50% longer
#define MAX_FREQ (1<<16)-32

template <int NSYM>
struct SIMPLE_MODEL {
    enum { STEP=8 };

    SIMPLE_MODEL();
    inline void encodeSymbol(RangeCoder *rc, uint16_t sym);
    inline int encodeNearSymbol(RangeCoder *rc, uint16_t sym, int dist);
    inline uint16_t decodeSymbol(RangeCoder *rc);

protected:
    void   normalize();

    uint32_t TotFreq;  // Total frequency
    uint32_t BubCnt;   // Periodic counter for bubble sort step

    // Array of Symbols approximately sorted by Freq. 
	struct SymFreqs {
		uint16_t Symbol;
		uint16_t Freq;
    } sentinel, F[NSYM+1];
};


template <int NSYM>
SIMPLE_MODEL<NSYM>::SIMPLE_MODEL() {
    for ( int i=0; i<NSYM; i++ ) {
		F[i].Symbol = i;
		F[i].Freq   = 1;
    }

    TotFreq         = NSYM;
    sentinel.Symbol = 0;
    sentinel.Freq   = MAX_FREQ; // Always first; simplifies sorting.
    BubCnt          = 0;

    F[NSYM].Freq = 0; // terminates normalize() loop. See below.
}


template <int NSYM>
void SIMPLE_MODEL<NSYM>::normalize() {
    /* Faster than F[i].Freq for 0 <= i < NSYM */
    TotFreq=0;
    for (SymFreqs *s = F; s->Freq; s++) {
		s->Freq -= s->Freq>>1;
		TotFreq += s->Freq;
    }
}

template <int NSYM>
inline void SIMPLE_MODEL<NSYM>::encodeSymbol(RangeCoder *rc, uint16_t sym) {
    SymFreqs *s = F;
    uint32_t AccFreq  = 0;

    while (s->Symbol != sym)
    	AccFreq += s++->Freq;

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
    	normalize();

    /* Keep approx sorted */
    if (((++BubCnt & 15) == 0) && s[0].Freq > s[-1].Freq) {
		SymFreqs t = s[0];
		s[0]  = s[-1];
		s[-1] = t;
    }
}

template <int NSYM>
inline int SIMPLE_MODEL<NSYM>::encodeNearSymbol(RangeCoder *rc, uint16_t sym, int dist) {
    SymFreqs *s = F;
    uint32_t AccFreq  = 0;
    int new_sym;

    while ( ABS((signed int)s->Symbol - (signed int)sym) > dist )
    	AccFreq += s++->Freq;

    rc->Encode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
    	normalize();

    new_sym = s->Symbol;

    /* Keep approx sorted */
    if (((++BubCnt&15)==0) && s[0].Freq > s[-1].Freq) {
		SymFreqs t = s[0];
		s[0] = s[-1];
		s[-1] = t;
    }

    return new_sym;
}

template <int NSYM>
inline uint16_t SIMPLE_MODEL<NSYM>::decodeSymbol(RangeCoder *rc) {
    SymFreqs* s = F;
    uint32_t freq = rc->GetFreq(TotFreq);
    uint32_t AccFreq;

    for (AccFreq = 0; (AccFreq += s->Freq) <= freq; s++);
    	AccFreq -= s->Freq;

    rc->Decode(AccFreq, s->Freq, TotFreq);
    s->Freq += STEP;
    TotFreq += STEP;

    if (TotFreq > MAX_FREQ)
    	normalize();

    /* Keep approx sorted */
    if (((++BubCnt&15)==0) && s[0].Freq > s[-1].Freq) {
		SymFreqs t = s[0];
		s[0] = s[-1];
		s[-1] = t;
		return t.Symbol;
    }

    return s->Symbol;
}

// Start main

/* QBITS is the size of the quality context, 2x quals */
#define QBITS 12
#define QSIZE (1<<QBITS)

/* Keep as a power of 2 */
//#define QMAX 128
#define QMAX 64

/*
 * SSE support to allow use of memory prefetching. It's only minor, but
 * it all helps.
 *
 * With    on -s6 -q3 -b (40million lines)
 *    encode: 2m55.691s+0m8.813s
 *    decode: 3m19.492s+0m2.720s
 *
 * Without on -s6 -q3 -b
 *    encode: 2m57.115s+0m3.260s
 *    decode: 3m46.338s+0m2.856s
 */
#ifdef __SSE__
#   include <xmmintrin.h>
#else
#   define _mm_prefetch(a,b)
#endif

struct fqz {
	fqz() : qlevel(2), model_qual(nullptr){
		int qsize = QSIZE;
		if (qlevel > 1) qsize *= 16;
		if (qlevel > 2) qsize *= 16;

		model_qual = new SIMPLE_MODEL<QMAX>[qsize];
	}

	~fqz(){ delete [] model_qual; }

    // Quality
    SIMPLE_MODEL<QMAX>* model_qual;
#define SMALL_QMASK (QSIZE-1)

    void encode_qual(RangeCoder *rc, char *qual, int len);
    void decode_qual(RangeCoder *rc, char *qual, int len);

    int qlevel;

};


