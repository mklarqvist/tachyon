#ifndef CORE_BASE_ENTRYHOTMETA_H_
#define CORE_BASE_ENTRYHOTMETA_H_

#include "../../io/BasicBuffer.h"

namespace Tachyon{
namespace Core{
/*
 TomahawkEntryMetaBase is used for reinterpreting
 byte streams of TomahawkEntryMeta
 Number of runs can be inferred from the sample
 number and byte length of the stream
 */
struct __attribute__((packed)) EntryHotMeta{
	typedef EntryHotMeta self_type;
	typedef IO::BasicBuffer buffer_type;

public:
	typedef struct __meta_controller{
		explicit __meta_controller(void) :
				anyMissing(0),
				allPhased(0),
				mixed_phasing(0),
				biallelic(0),
				simple(0),
				rle(0),
				rle_type(0),
				diploid(0),
				mixed_ploidy(0),
				unused(6)
		{}
		~__meta_controller(){}

		U16  anyMissing: 1,  // any missing
             allPhased: 1,   // all phased
			 mixed_phasing: 1,// has mixed phasing
			 biallelic: 1,   // is biallelic
			 simple: 1,      // is simple SNV->SNV
			 rle: 1,         // uses RLE compression
			 rle_type: 2,   // type of RLE (BYTE, U16, U32, U64)
			 diploid: 1,    // is diploid
			 mixed_ploidy: 1, // has mixed ploidy (e.g. X chromosome or CNV)
			 unused: 6; // reserved
	} controller_type;

public:
	EntryHotMeta() :
		position(0),
		ref_alt(0),
		AF(0),
		FILTER_map_ID(0),
		INFO_map_ID(0),
		FORMAT_map_ID(0),
		virtual_offset_cold_meta(0),
		virtual_offset_gt(0),
		n_runs(0)
	{}
	~EntryHotMeta(){}

	inline const bool isSingleton(void) const{ return(this->AF == 0); }
	inline const bool isSimpleSNV(void) const{ return(this->controller.biallelic == true && this->controller.simple == true); }
	inline const bool isRLE(void) const{ return(this->controller.rle); }
	inline const bool isDiploid(void) const{ return(this->controller.diploid); }
	inline const bool isMixedPloidy(void) const{ return(this->controller.mixed_ploidy); }

	inline const U32& getRuns(void) const{ return(this->n_runs); }
	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << entry.position << '\t' <<
			   (int)*reinterpret_cast<const BYTE* const>(&entry.controller) << '\t' <<
			   (int)entry.ref_alt << '\t' <<
			   entry.AF << '\t' <<
			   entry.virtual_offset_cold_meta << '\t' <<
			   entry.virtual_offset_gt << '\t' <<
			   entry.n_runs;
		return(out);
	}

	// Overload operator+= for basic buffer
	friend buffer_type& operator+=(buffer_type& buffer, const self_type& entry){
		buffer += *reinterpret_cast<const U16* const>(&entry.controller);
		buffer += entry.position;
		buffer += entry.ref_alt;
		buffer += entry.AF;
		buffer += entry.FILTER_map_ID;
		buffer += entry.INFO_map_ID;
		buffer += entry.FORMAT_map_ID;
		buffer += entry.virtual_offset_cold_meta;
		buffer += entry.virtual_offset_gt;
		buffer += entry.n_runs;
		return(buffer);
	}

public:
	controller_type controller;
	U32 position; // is block.minPosition + position

	// Most sites are bi-allelic and simple SNVs
	// sites that are not bi-allelic and not simple
	// will be encoded in the complex meta section
	BYTE ref_alt;

	// AF is pre-computed
	float AF;

	// FILTER map
	U16 FILTER_map_ID;
	// INFO map
	U16 INFO_map_ID;
	// FORMAT map
	U16 FORMAT_map_ID;

	// Hot-cold split structure. pointer to cold data
	// since a pointer cannot be read from a byte
	// stream as its memory location changes
	// we have to provide the pointer as an ABSOLUTE
	// virtual stream offset relative to the complex
	// start position into the complex byte stream
	U32 virtual_offset_cold_meta;
	// offset to the byte end of this entry in the stream
	// (either RLE or simple stream depending on context)
	// this allows fast iteration when switching between
	// the two compression approaches
	U32 virtual_offset_gt;
	// Runs
	U32 n_runs;
};

}
}

#endif /* CORE_BASE_ENTRYHOTMETA_H_ */
