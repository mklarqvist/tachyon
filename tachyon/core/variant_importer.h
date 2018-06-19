#ifndef CORE_VARIANT_IMPORTER_H_
#define CORE_VARIANT_IMPORTER_H_

#include "../algorithm/compression/compression_manager.h"
#include "../algorithm/compression/genotype_encoder.h"
#include "../algorithm/permutation/radix_sort_gt.h"
#include "../support/type_definitions.h"
#include "../support/helpers.h"
#include "../io/bcf/BCFReader.h"
#include "../containers/variant_block.h"
#include "../index/index_entry.h"
#include "../index/index_index_entry.h"
#include "variant_importer_container_stats.h"
#include "../algorithm/timer.h"
#include "variant_import_writer.h"

namespace tachyon {

/**<
 * Settings for VariantImporter class
 */
struct VariantImporterSettings{
public:
	VariantImporterSettings() :
		permute_genotypes(true),
		encrypt_data(false),
		drop_invariant_sites(false),
		checkpoint_n_snps(1000),
		checkpoint_bases(10e6),
		n_threads(std::thread::hardware_concurrency()),
		info_end_key(-1),
		info_svlen_key(-1),
		compression_level(6)
	{}

	~VariantImporterSettings() = default;

	std::string getInterpretedString(void) const{
		return(std::string("##tachyon_importInterpretedCommand=input_file=" + this->input_file +
		   ";output_prefix=" + this->output_prefix +
		   ";checkpoint_snps=" + std::to_string(this->checkpoint_n_snps) +
		   ";checkpoint_bases=" + std::to_string(this->checkpoint_bases) +
		   ";compression_level=" + std::to_string(this->compression_level)
		));
	}

	inline void setInputFile(const std::string& input_name){ this->input_file = input_name; }
	inline void setOutputPrefix(const std::string& output_prefix){ this->output_prefix = output_prefix; }
	inline void setThreads(const U32 n_threads){ this->n_threads = n_threads; }
	inline void setPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void setEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void setCompressionLevel(const U32 compression_level){ this->compression_level = compression_level; }

public:
	bool permute_genotypes;   // permute GT flag
	bool encrypt_data;        // encryption flag
	bool drop_invariant_sites;// drop sites that are invariant
	U32 checkpoint_n_snps;    // number of variants until checkpointing
	U32 checkpoint_bases;     // number of bases until checkpointing
	U32 n_threads;            // number of parallel importer threads
	S32 info_end_key;         // key mapping to the INFO field END
	S32 info_svlen_key;       // key mapping to the INFO field SVLEN
	U32 compression_level;    // compression level sent to ZSTD
	std::string input_file;   // input file name
	std::string output_prefix;// output file prefix
};

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
	typedef VariantImporterSettings         settings_type;

public:
	VariantImporter();
	VariantImporter(const settings_type& settings);

	VariantImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint_size, const double checkpoint_bases);
	~VariantImporter();
	bool Build();

	void setWriterTypeFile(void){ this->writer = new writer_file_type; }
	void setWriterTypeStream(void){ this->writer = new writer_stream_type; }

private:
	bool BuildBCF();  // import a BCF file
	bool addSite(meta_type& meta, bcf_entry_type& line); // Import a BCF line
	bool addGenotypes(bcf_reader_type& bcf_reader, meta_type* meta_entries);
	bool parseBCFBody(meta_type& meta, bcf_entry_type& line);

private:
	settings_type settings_; // internal settings
	bool GT_available_;

	// Stats
	import_stats_type stats_basic;
	import_stats_type stats_info;
	import_stats_type stats_format;

	// Read/write fields
	writer_interface_type* writer;      // writer

	index_entry_type  index_entry; // streaming index entry
	radix_sorter_type permutator;  // GT permuter
	header_type*      header;      // header
	gt_encoder_type   encoder;     // RLE packer

	compression_manager_type compression_manager;

	// Data container
	block_type block;
};


} /* namespace Tachyon */

#endif /* CORE_VARIANT_IMPORTER_H_ */
