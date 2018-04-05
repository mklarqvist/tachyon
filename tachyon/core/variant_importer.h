#ifndef CORE_VARIANT_IMPORTER_H_
#define CORE_VARIANT_IMPORTER_H_

#include "../algorithm/compression/compression_manager.h"
#include "../algorithm/compression/genotype_encoder.h"
#include "../algorithm/permutation/radix_sort_gt.h"
#include "../support/type_definitions.h"
#include "../support/helpers.h"
#include "../io/bcf/BCFReader.h"
#include "../containers/variantblock.h"
#include "../index/index_entry.h"
#include "../index/index_index_entry.h"
#include "variant_importer_container_stats.h"
#include "../algorithm/timer.h"
#include "variant_import_writer.h"

namespace tachyon {

class VariantImporter {
private:
	typedef VariantImporter                 self_type;
	typedef VariantImportWriterInterface    writer_interface_type;
	typedef VariantImportWriterFile         writer_file_type;
	typedef VariantImportWriterStream       writer_stream_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef vcf::VCFHeader                  header_type;
	typedef index::IndexEntry               index_entry_type;
	typedef bcf::BCFReader                  bcf_reader_type;
	typedef bcf::BCFEntry                   bcf_entry_type;
	typedef algorithm::CompressionManager   compression_manager_type;
	typedef algorithm::RadixSortGT          radix_sorter_type;
	typedef algorithm::PermutationManager   permutation_type;
	typedef algorithm::GenotypeEncoder      gt_encoder_type;
	typedef containers::DataContainer       stream_container;
	typedef containers::HashContainer       hash_container_type;
	typedef containers::HashVectorContainer hash_vector_container_type;
	typedef containers::VariantBlock        block_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef core::MetaEntry                 meta_type;

public:
	VariantImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint_size, const double checkpoint_bases);
	~VariantImporter();
	bool Build();

	inline void setPermute(const bool yes){ this->permute = yes; }
	inline void setEncrypt(const bool yes){ this->encrypt = yes; }
	void setWriterTypeFile(void){ this->writer = new writer_file_type; }
	void setWriterTypeStream(void){ this->writer = new writer_stream_type; }

private:
	bool BuildBCF();  // import a BCF file
	bool add(bcf_entry_type& line); // Import a BCF line
	bool parseBCFBody(meta_type& meta, bcf_entry_type& line);

private:
	bool GT_available_;
	bool permute;            // permute GT flag
	bool encrypt;            // encryption flag
	U32 checkpoint_n_snps;   // number of variants until checkpointing
	double checkpoint_bases; // number of bases until checkpointing

	// Stats
	import_stats_type stats_basic;
	import_stats_type stats_info;
	import_stats_type stats_format;

	// Read/write fields
	std::string inputFile;   // input file name
	std::string outputPrefix;// output file prefix
	writer_interface_type* writer;      // writer

	index_entry_type  index_entry; // streaming index entry
	radix_sorter_type permutator;  // GT permuter
	header_type*      header;      // header
	gt_encoder_type   encoder;     // RLE packer

	compression_manager_type compression_manager;

	// Data container
	block_type block;

	// temp
	//algorithm::GenotypeNearestNeighbour* nn;

};


} /* namespace Tachyon */

#endif /* CORE_VARIANT_IMPORTER_H_ */
