#ifndef CORE_IMPORTER_H_
#define CORE_IMPORTER_H_

#include "../algorithm/compression/compression_manager.h"
#include "../algorithm/compression/genotype_encoder.h"
#include "../algorithm/permutation/radix_sort_gt.h"
#include "../support/type_definitions.h"
#include "../support/helpers.h"
#include "../io/bcf/BCFReader.h"
#include "ImportWriter.h"
#include "ImporterStats.h"
#include "../containers/DataBlock.h"
#include "../index/index_entry.h"
#include "../index/index_index_entry.h"

namespace tachyon {

class Importer {
private:
	typedef Importer                        self_type;
	typedef io::BasicReader                 reader_type;
	typedef vcf::VCFHeader                  header_type;
	typedef ImportWriter                    writer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef algorithm::GenotypeEncoder      gt_encoder_type;
	typedef index::IndexEntry               index_entry_type;
	typedef bcf::BCFReader                  bcf_reader_type;
	typedef bcf::BCFEntry                   bcf_entry_type;
	typedef algorithm::RadixSortGT          radix_sorter_type;
	typedef containers::DataContainer       stream_container;
	typedef algorithm::PermutationManager   permutation_type;
	typedef core::MetaHot                   meta_type;
	typedef containers::HashContainer       hash_container_type;
	typedef containers::HashVectorContainer hash_vector_container_type;
	typedef containers::DataBlock           block_type;
	typedef support::ImporterStats          import_stats_type;
	typedef algorithm::CompressionManager   compression_manager_type;

public:
	Importer(std::string inputFile, std::string outputPrefix, const U32 checkpoint_size, const double checkpoint_bases);
	~Importer();
	bool Build();

	inline void setPermute(const bool yes){ this->permute = yes; }

private:
	bool BuildBCF();  // import a BCF file
	bool parseBCFLine(bcf_entry_type& line); // Import a BCF line
	bool parseBCFBody(meta_type& meta, bcf_entry_type& line);
	void resetHashes(void);

private:
	bool permute;               // permute GT flag
	U32 checkpoint_n_snps;      // number of variants until checkpointing
	double checkpoint_bases;    // number of bases until checkpointing

	// Stats
	import_stats_type import_uncompressed_stats;
	import_stats_type import_compressed_stats;

	// Read/write fields
	std::string inputFile;   // input file name
	std::string outputPrefix;// output file prefix
	reader_type reader;      // reader
	writer_type writer;      // writer

	index_entry_type  index_entry;  // Header index
	radix_sorter_type permutator;
	header_type*      header;     // header
	gt_encoder_type   encoder;     // RLE packer

	compression_manager_type compression_manager;

	// Data container
	block_type block;

	// Use during import only
	hash_container_type info_fields;
	hash_container_type format_fields;
	hash_container_type filter_fields;
	hash_vector_container_type info_patterns;
	hash_vector_container_type format_patterns;
	hash_vector_container_type filter_patterns;
};


} /* namespace Tachyon */

#endif /* CORE_IMPORTER_H_ */
