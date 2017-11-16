#ifndef TOMAHAWKIMPORTER_H_
#define TOMAHAWKIMPORTER_H_

#include "../algorithm/compression/EncoderGenotypesRLE.h"
#include "../index/IndexEntry.h"
#include "../index/IndexBlockEntry.h"
#include "../algorithm/permutation/RadixSortGT.h"
#include "TomahawkImportWriter.h"

namespace Tomahawk {

class TomahawkImporter {
	typedef TomahawkImporter self_type;
	typedef reader reader_type;
	typedef VCF::VCFHeader header_type;
	typedef TomahawkImportWriter writer_type;
	typedef VCF::VCFLine line_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Algorithm::EncoderGenotypesRLE encoder_type;
	typedef Totempole::IndexEntry totempole_entry_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry bcf_entry_type;
	typedef Algorithm::RadixSortGT radix_sorter_type;

	/*
	 This supportive structure keeps track of the current and
	 previous contig identifiers and the previous obseved position.
	 This information is necessary to guarantee the sort-order of
	 the output Tomahawk file required for indexing.
	 Note that contigID is a pointer as this is required by our
	 hash-table implementation as a return value
	 */
	struct __InternalHelper{
		__InternalHelper():
			contigID(nullptr),
			prevcontigID(-1),
			previous_position(-1)
		{}
		S32* contigID;			// current contigID
		S32 prevcontigID;		// previous contigID
		S32 previous_position;	// current position
	} sort_order_helper;

public:
	TomahawkImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint_size);
	~TomahawkImporter();
	bool Build();

private:
	bool BuildBCF();  // import a BCF file
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line

private:
	bool permutateData(bcf_reader_type& reader);

private:
	U32 checkpoint_size;      // number of variants until checkpointing
	U32 block_flush_limit;    // limit in bytes when to flush to disk
	std::string inputFile;    // input file name
	std::string outputPrefix; // output file prefix
	reader_type reader_;      // reader
	writer_type writer_;      // writer
	buffer_type meta_buffer;  // meta buffer
	buffer_type encode_rle_buffer;   // RLE buffer
	buffer_type encode_simple_buffer;   // RLE buffer
	totempole_entry_type totempole_entry;  // totempole entry for indexing
	radix_sorter_type permutator;
	header_type* header_;     // header
	encoder_type* encoder;   // RLE packer
};


} /* namespace Tomahawk */

#endif /* TOMAHAWKIMPORTER_H_ */
