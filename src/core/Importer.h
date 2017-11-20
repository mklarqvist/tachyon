#ifndef CORE_IMPORTER_H_
#define CORE_IMPORTER_H_

#include "../support/TypeDefinitions.h"
#include "../support/helpers.h"
#include "../io/bcf/BCFReader.h"
#include "BlockEntry.h"
#include "ImportWriter.h"
#include "../algorithm/permutation/RadixSortGT.h"

namespace Tachyon {

class Importer {
	typedef Importer self_type;
	typedef reader reader_type;
	typedef VCF::VCFHeader header_type;
	typedef ImportWriter writer_type;
	typedef IO::BasicBuffer buffer_type;
	typedef Algorithm::EncoderGenotypesRLE encoder_type;
	typedef Totempole::IndexEntry totempole_entry_type;
	typedef BCF::BCFReader bcf_reader_type;
	typedef BCF::BCFEntry bcf_entry_type;
	typedef Algorithm::RadixSortGT radix_sorter_type;
	typedef Core::StreamContainer stream_container;
	typedef Core::PermutationManager permutation_type;
	typedef Core::EntryHotMetaBase meta_type;
	typedef Core::Support::HashContainer hash_container_type;
	typedef Core::Support::HashVectorContainer hash_vector_container_type;
	typedef Core::BlockEntry block_type;

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
	Importer(std::string inputFile, std::string outputPrefix, const U32 checkpoint_size);
	~Importer();
	bool Build();

private:
	bool BuildBCF();  // import a BCF file
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line
	bool parseBCFBody(meta_type& meta, bcf_entry_type& line);

private:
	bool permutateData(bcf_reader_type& reader);
	void resetHashes(void);

private:
	U32 checkpoint_size;      // number of variants until checkpointing
	std::string inputFile;    // input file name
	std::string outputPrefix; // output file prefix
	reader_type reader_;      // reader
	writer_type writer_;      // writer

	totempole_entry_type totempole_entry;  // totempole entry for indexing
	radix_sorter_type permutator;
	header_type* header_;     // header
	encoder_type encoder;     // RLE packer

	block_type block;

	hash_container_type info_fields;
	hash_container_type format_fields;
	hash_container_type filter_fields;
	hash_vector_container_type info_patterns;
	hash_vector_container_type format_patterns;
	hash_vector_container_type filter_patterns;
};


} /* namespace Tachyon */

#endif /* CORE_IMPORTER_H_ */
