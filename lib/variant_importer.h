#ifndef VARIANT_IMPORTER_H_
#define VARIANT_IMPORTER_H_

#include <unordered_map>

#include "algorithm/compression/compression_manager.h"
#include "algorithm/compression/genotype_encoder.h"
#include "algorithm/timer.h"
#include "containers/variant_block.h"
#include "io/variant_import_writer.h"
#include "core/variant_importer_container_stats.h"
#include "index/variant_index_entry.h"
#include "index/variant_index_meta_entry.h"
#include "io/vcf_utils.h"
#include "support/helpers.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "core/footer/footer.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "algorithm/permutation/genotype_sorter.h"

#include "core/variant_import_settings.h"

namespace tachyon {

class VariantImporter {
public:
	typedef VariantImporter                 self_type;
	typedef VariantWriterInterface    writer_interface_type;
	typedef VariantWriterFile         writer_file_type;
	typedef VariantWriterStream       writer_stream_type;
	typedef io::VcfReader                   vcf_reader_type;
	typedef algorithm::CompressionManager   compression_manager_type;
	typedef VariantImporterSettings         settings_type;
	typedef std::unordered_map<uint32_t, uint32_t>    reorder_map_type;
	typedef std::unordered_map<uint64_t, uint32_t>    hash_map_type;

public:
	VariantImporter();
	VariantImporter(const settings_type& settings);
	~VariantImporter();

	bool Build();
	bool BuildParallel();

	void clear(void);

	inline void SetWriterTypeFile(void)  { this->writer = new writer_file_type;   }
	inline void SetWriterTypeStream(void){ this->writer = new writer_stream_type; }

private:
	bool WriteFinal(algorithm::VariantDigestManager& checksums);
	bool WriteKeychain(const encryption::Keychain<>& keychain);
	bool WriteYonHeader();
	void UpdateHeaderImport(VariantHeader& header);

private:
	settings_type settings_; // import settings

	// Read/write fields
	writer_interface_type* writer; // writer

	compression_manager_type compression_manager;

	// Map from BCF global FORMAT/INFO/FILTER IDX to local IDX such that
	// FORMAT maps to [0, f-1], and INFO maps to [0, i-1] and FILTER to
	// [0,l-1] and where f+i+l = n, where n is the total number of fields.
	//
	//                    Global    Local
	// std::unordered_map<uint32_t, uint32_t> filter_reorder_map_;
	reorder_map_type filter_reorder_map_;
	reorder_map_type info_reorder_map_;
	reorder_map_type format_reorder_map_;
	reorder_map_type contig_reorder_map_;

	std::unique_ptr<vcf_reader_type> vcf_reader_;
	VariantHeader yon_header_;

	hash_map_type block_hash_map;
};

} /* namespace Tachyon */

#endif /* VARIANT_IMPORTER_H_ */
