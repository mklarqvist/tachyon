#ifndef VARIANT_IMPORTER_H_
#define VARIANT_IMPORTER_H_

#include <unordered_map>
#include <mutex>
#include <condition_variable>

#include "algorithm/compression/compression_manager.h"
#include "algorithm/compression/genotype_encoder.h"
#include "algorithm/timer.h"
#include "containers/variant_block.h"
#include "core/variant_import_writer.h"
#include "core/variant_importer_container_stats.h"
#include "index/variant_index_entry.h"
#include "index/variant_index_meta_entry.h"
#include "io/vcf_utils.h"
#include "support/helpers.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "core/footer/footer.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "algorithm/permutation/genotype_sorter.h"

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
		compression_level(6),
		htslib_extra_threads(0)
	{}

	~VariantImporterSettings() = default;

	std::string GetInterpretedString(void) const{
		return(std::string("{\"input_file\":\"" + this->input_file +
		   "\",\"output_prefix\":\"" + this->output_prefix +
		   "\",\"checkpoint_snps\":" + std::to_string(this->checkpoint_n_snps) +
		   ",\"checkpoint_bases\":" + std::to_string(this->checkpoint_bases) +
		   ",\"compression_level\":" + std::to_string(this->compression_level) +
		   "}"
		));
	}

	inline void SetInputFile(const std::string& input_name){ this->input_file = input_name; }
	inline void SetOutputPrefix(const std::string& output_prefix){ this->output_prefix = output_prefix; }
	inline void SetThreads(const uint32_t n_threads){ this->n_threads = n_threads; }
	inline void SetPermute(const bool yes){ this->permute_genotypes = yes; }
	inline void SetEncrypt(const bool yes){ this->encrypt_data = yes; }
	inline void SetCompressionLevel(const uint32_t compression_level){ this->compression_level = compression_level; }

public:
	bool permute_genotypes;   // permute GT flag
	bool encrypt_data;        // encryption flag
	bool drop_invariant_sites;// drop sites that are invariant
	uint32_t checkpoint_n_snps;    // number of variants until checkpointing
	uint32_t checkpoint_bases;     // number of bases until checkpointing
	uint32_t n_threads;            // number of parallel importer threads
	int32_t info_end_key;         // key mapping to the INFO field END
	int32_t info_svlen_key;       // key mapping to the INFO field SVLEN
	uint32_t compression_level;    // compression level sent to ZSTD
	std::string input_file;   // input file name
	std::string output_prefix;// output file prefix
	uint32_t htslib_extra_threads;
};

/**<
 * Cyclic queue payload object used when importing htslib Vcf
 * entries. Used exclusively in the yon_pool_vcf struct.
 */
struct yon_pool_vcfc_payload {
	yon_pool_vcfc_payload() : block_id(0), c(nullptr){}
	yon_pool_vcfc_payload(const uint32_t bid, containers::VcfContainer* vc) : block_id(bid), c(vc){}
	~yon_pool_vcfc_payload(){ delete c; }
	yon_pool_vcfc_payload(const yon_pool_vcfc_payload& other) = delete; // copy is not allowed
	yon_pool_vcfc_payload& operator=(const yon_pool_vcfc_payload& other) = delete; // copy assign is not allowed
	yon_pool_vcfc_payload(yon_pool_vcfc_payload&& other) noexcept : block_id(other.block_id), c(nullptr){
		std::swap(this->c, other.c);
	}
	yon_pool_vcfc_payload& operator=(yon_pool_vcfc_payload&& other) noexcept{
		if(this == &other) return(*this);

		this->block_id = other.block_id;
		delete this->c;
		this->c = nullptr;
		std::swap(this->c, other.c);
		return(*this);
	}

	uint32_t block_id;
	containers::VcfContainer* c;
};

/**<
 * Cyclic pool of loaded VcfContainers. Uses conditional variables
 * and mutex locks to control flow to and from this object. Use
 * emplace() to add a payload and pop() to retrieve an object.
 */
struct yon_pool_vcfc {
public:
	yon_pool_vcfc(void) :
		c(nullptr),
		n_capacity(0), n_c(0),
		front(0), rear(0), alive(false)
	{

	}

	yon_pool_vcfc(uint32_t capacity) :
		c(new yon_pool_vcfc_payload*[capacity]),
		n_capacity(capacity), n_c(0),
		front(0), rear(0), alive(false)
	{

	}

	yon_pool_vcfc(const yon_pool_vcfc& other) = delete; // disallow copy
	yon_pool_vcfc& operator=(const yon_pool_vcfc& other) = delete; // disallow assign copy

	yon_pool_vcfc(yon_pool_vcfc&& other) noexcept :
		c(nullptr),
		n_capacity(other.n_capacity), n_c(other.n_c),
		front(other.front), rear(other.rear), alive(other.alive.load())
	{
		std::swap(this->c, other.c);
	}

	yon_pool_vcfc& operator=(yon_pool_vcfc&& other) noexcept {
		if(this == &other){ return(*this); }
		delete [] this->c;
		this->c = nullptr;
		this->n_capacity = other.n_capacity;
		this->n_c = other.n_c;
		this->front = other.front;
		this->rear = other.rear;
		this->alive = other.alive.load();
		std::swap(this->c, other.c);
		return(*this);
	}

	~yon_pool_vcfc(){
		delete [] this->c;
	}

	/**<
	 * Add payload to the cyclic queue. If the queue is full then
	 * wait until an item has been popped.
	 * @param data Input pointer to payload.
	 */
	void emplace(yon_pool_vcfc_payload* data){
		std::unique_lock<std::mutex> l(lock);

		this->not_full.wait(l, [this](){
			if(this->n_c == 0 && this->alive == false) return true;
			return this->n_c != this->n_capacity;
		});

		if(this->n_c == 0 && this->alive == false){
			//std::cerr << "is empty in wait full" << std::endl;
			l.unlock();
			this->not_full.notify_all();
			this->not_empty.notify_all();
			return;
		}

		this->c[this->rear] = data;
		this->rear = (this->rear + 1) % this->n_capacity;
		++this->n_c;

		l.unlock();
		this->not_empty.notify_one();
	}

	/**<
	 * Retrieve payload from the cyclic queue. If the queue is empty
	 * then wait until an item has been inserted.
	 * @return Returns a pointer to the retrieved payload or a nullpointer in the special case the producer operations has finished.
	 */
	yon_pool_vcfc_payload* pop(void){
		std::unique_lock<std::mutex> l(lock);

		this->not_empty.wait(l, [this](){
			if(this->n_c == 0 && this->alive == false) return true;
			return this->n_c != 0;
		});

		if(this->n_c == 0){
			//std::cerr << "is empty return nullptr" << std::endl;
			l.unlock();
			this->not_empty.notify_all();
			this->not_full.notify_all();
			return nullptr;
		}

		yon_pool_vcfc_payload* result = this->c[this->front];
		this->c[this->front] = nullptr;
		this->front = (this->front + 1) % this->n_capacity;
		--this->n_c;

		l.unlock();
		this->not_full.notify_one();

		return result;
	}

public:
	yon_pool_vcfc_payload** c;
	int n_capacity;
	int n_c;
	int front;
	int rear;
	std::atomic<bool> alive;
	std::mutex lock;
	std::condition_variable not_full;
	std::condition_variable not_empty;
};

/**<
 * Single producer (albeit spawning many htslib decompressor threads)
 * reading htslib Vcf records from a stream into a VcfContainer and
 * inserts that object into the cyclic producer queue.
 */
struct yon_producer_vcfc {
public:
	yon_producer_vcfc(const VariantImporterSettings& settings,
	                  std::unique_ptr<io::VcfReader>& reader,
	                  uint32_t pool_size) :
		all_finished(false),
		n_rcds_loaded(0),
		data_available(false),
		data_pool(pool_size),
		settings(settings),
		reader(reader)
	{}

	/**<
	 * Spawn worker reading data into the producer data pool.
	 * @return Returns a reference to the spawned thread.
	 */
	std::thread& Start(void){
		this->thread_ = std::thread(&yon_producer_vcfc::Produce, this);
		return(this->thread_);
	}

private:
	bool Produce(void){
		uint32_t n_blocks = 0;
		this->data_available  = true;
		this->data_pool.alive = true;
		while(true){
			// Retrieve bcf1_t records using htslib and lazy evaluate them. Stop
			// after retrieving a set number of variants or if the interval between
			// the smallest and largest variant exceeds some distance in base pairs.
			// Vcf records are NOT lazy evaluated (unpacked) in the producer step
			// that is the job of the consumers.
			if(this->container.GetVariants(this->settings.checkpoint_n_snps,
			                               this->settings.checkpoint_bases,
			                               this->reader,
			                               0) == false)
			{

				this->data_available  = false;
				this->data_pool.alive = false;
				while(data_pool.n_c && all_finished == false){
					data_pool.not_empty.notify_all();
					data_pool.not_full.notify_all();
				}

				break;
			}
			this->n_rcds_loaded += this->container.sizeWithoutCarryOver();
			// Peculiar syntax for adding a new payload to the cyclic queue. The
			// move semantics is required for intended functionality!
			this->data_pool.emplace(new yon_pool_vcfc_payload(n_blocks++, new containers::VcfContainer(std::move(this->container))));
		}
		return true;
	}

public:
	bool all_finished;
	uint64_t n_rcds_loaded;
	std::atomic<bool> data_available;
	yon_pool_vcfc data_pool;
	containers::VcfContainer container;
	const VariantImporterSettings& settings;
	std::thread thread_;
	std::unique_ptr<io::VcfReader>& reader;
};

/**<
 * Slave worker instantiation of VariantImporter. Used exclusively during htslib
 * Vcf importing. Data has to be provided through the cyclic queue generated by the
 * appropriate producer.
 */
class VcfImporter {
public:
	typedef VcfImporter                     self_type;
	typedef VariantImportWriterInterface    writer_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef index::VariantIndexEntry        index_entry_type;
	typedef io::VcfHeader                   vcf_header_type;
	typedef containers::VcfContainer        vcf_container_type;
	typedef algorithm::CompressionManager   compression_manager_type;
	typedef algorithm::GenotypeSorter       radix_sorter_type;
	typedef algorithm::GenotypeEncoder      gt_encoder_type;
	typedef containers::DataContainer       stream_container;
	typedef containers::VariantBlock        block_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef core::MetaEntry                 meta_type;
	typedef VariantImporterSettings         settings_type;
	typedef std::unordered_map<uint32_t, uint32_t> reorder_map_type;
	typedef std::unordered_map<uint64_t, uint32_t> hash_map_type;

public:
	VcfImporter() : n_blocks_processed(0), block_id(0), vcf_header_(nullptr), GT_available_(false){}
	VcfImporter(const settings_type& settings) : n_blocks_processed(0), block_id(0), vcf_header_(nullptr), GT_available_(false){}
	VcfImporter(const self_type& other) = delete;
	VcfImporter(self_type&& other) noexcept = delete;
	VcfImporter& operator=(const self_type& other) = delete;
	VcfImporter& operator=(self_type&& other) = delete;
	~VcfImporter(){
		// do not delete vcf_header, it is not owned by this class
	}

	VcfImporter& operator+=(const VcfImporter& other){
		this->n_blocks_processed += other.n_blocks_processed;
		this->index += other.index;
		return(*this);
	}

	/**<
	 * Import the contents from a provided VcfContainer and target block id
	 * into the local VariantBlock while preparing it for writing.
	 * @param container
	 * @param block_id
	 * @return
	 */
	bool Add(vcf_container_type& container, const uint32_t block_id){
		// Clear current data.
		this->block.clear();
		this->index_entry.reset();

		// Assign new block id. This is provided by the producer.
		this->block_id = block_id;
		this->index_entry.blockID = block_id;

		if(this->n_blocks_processed++ == 0){
			// Allocate containers and offsets for this file.
			// This is not strictly necessary but prevents nasty resize
			// calls in most cases.
			this->block.Allocate(this->vcf_header_->info_fields_.size(),
								 this->vcf_header_->format_fields_.size(),
								 this->vcf_header_->filter_fields_.size());

			// Resize containers
			const uint32_t resize_to = this->settings_.checkpoint_n_snps * sizeof(uint32_t) * 2; // small initial allocation
			this->block.resize(resize_to);
		}

		// This pointer here is borrowed from the PPA manager
		// during import stages. Do not destroy the target block
		// before finishing with this.
		this->block.gt_ppa = &this->permutator.permutation_array;

		if(this->GT_available_ && this->settings_.permute_genotypes){
			// Only store the permutation array if the number of samples
			// are greater then one (1).
			if(this->vcf_header_->GetNumberSamples() > 1){
				if(this->permutator.Build(container, *this->vcf_header_) == false)
					return false;

				this->block.header.controller.hasGTPermuted = true;
			}
		}

		if(this->AddRecords(container) == false) return false;

		this->block.header.controller.hasGT  = this->GT_available_;
		this->block.header.n_variants        = container.sizeWithoutCarryOver();
		this->block.UpdateContainers();

		// Perform compression using standard parameters.
		if(!this->compression_manager.Compress(this->block, this->settings_.compression_level, 6)){
			std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to compress..." << std::endl;
			return false;
		}

		// Encrypt the variant block if desired.
		if(this->settings_.encrypt_data){
			// Generate field-unique identifiers.
			//this->GenerateIdentifiers();

			// Start encryption.
			//this->block.header.controller.anyEncrypted = true;
			//if(!encryption_manager.encrypt(this->block, keychain, YON_ENCRYPTION_AES_256_GCM)){
			//	std::cerr << utility::timestamp("ERROR","COMPRESSION") << "Failed to encrypt..." << std::endl;
			//}
		}

		this->block.Finalize();
		// After all compression and writing is finished the header
		// offsets are themselves compressed and stored in the block.
		this->block.PackFooter(); // Pack footer into buffer.
		this->compression_manager.zstd_codec.Compress(block.footer_support);

		return true;
	}

	bool AddRecords(const vcf_container_type& container);
	bool AddRecord(const vcf_container_type& container, const uint32_t position, meta_type& meta);
	bool AddVcfInfo(const bcf1_t* record, meta_type& meta);
	bool AddVcfFormatInfo(const bcf1_t* record, meta_type& meta);
	bool AddVcfFilterInfo(const bcf1_t* record, meta_type& meta);
	bool IndexRecord(const bcf1_t* record, const meta_type& meta);
	bool AddVcfInfoPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddVcfFormatPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddVcfFilterPattern(const std::vector<int>& pattern, meta_type& meta);
	bool AddGenotypes(const vcf_container_type& container, meta_type* meta_entries);

	//bool GenerateIdentifiers(void);

	inline void SetVcfHeader(vcf_header_type* header){
		this->vcf_header_ = header;
		// Setup genotype permuter and genotype encoder.
		this->permutator.SetSamples(header->GetNumberSamples());
		this->encoder.SetSamples(header->GetNumberSamples());
	}

public:
	uint32_t n_blocks_processed;
	uint32_t block_id;
	index::Index index;
	vcf_header_type* vcf_header_;

	settings_type settings_; // internal settings
	bool GT_available_;

	// Read/write fields
	index_entry_type  index_entry; // streaming index entry
	radix_sorter_type permutator;  // GT permuter
	gt_encoder_type   encoder;     // RLE packer

	compression_manager_type compression_manager;

	// Data container
	block_type block;

	reorder_map_type filter_reorder_map_;
	reorder_map_type info_reorder_map_;
	reorder_map_type format_reorder_map_;
	reorder_map_type contig_reorder_map_;
};

class VariantImporter {
public:
	typedef VariantImporter                 self_type;
	typedef VariantImportWriterInterface    writer_interface_type;
	typedef VariantImportWriterFile         writer_file_type;
	typedef VariantImportWriterStream       writer_stream_type;
	typedef io::BasicBuffer                 buffer_type;
	typedef index::VariantIndexEntry        index_entry_type;
	typedef io::VcfReader                   vcf_reader_type;
	typedef containers::VcfContainer        vcf_container_type;
	typedef algorithm::CompressionManager   compression_manager_type;
	typedef algorithm::GenotypeSorter       radix_sorter_type;
	typedef algorithm::GenotypeEncoder      gt_encoder_type;
	typedef containers::DataContainer       stream_container;
	typedef containers::VariantBlock        block_type;
	typedef support::VariantImporterContainerStats import_stats_type;
	typedef core::MetaEntry                 meta_type;
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
	bool GT_available_;

	// Stats
	import_stats_type stats_basic;
	import_stats_type stats_info;
	import_stats_type stats_format;

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

	hash_map_type block_hash_map;
};

// Synchronised writer.
struct yon_writer_sync {
	yon_writer_sync() : n_written_rcds(0), next_block_id(0), alive(true), writer(nullptr){}
	~yon_writer_sync(){}

	void emplace(uint32_t block_id, const containers::VcfContainer& container, VcfImporter& importer){
		std::unique_lock<std::mutex> l(lock);

		cv_next_checkpoint.wait(l, [this, block_id](){
			if(this->alive == false) return true;
			return this->next_block_id == block_id;
		});

		if(this->alive == false){
			//std::cerr << "is exit condition in wblock" << std::endl;
			l.unlock();
			cv_next_checkpoint.notify_all();
			return;
		}

		// Write data container.
		this->Write(importer.block, container, importer.index_entry);
		// Update compression/storage statistics.
		importer.block.UpdateOutputStatistics(this->stats_basic, this->stats_info, this->stats_format);
		//std::cerr << utility::timestamp("LOG") << "Writing: " << this->next_block_id << ": " << importer.block.GetCompressedSize() << "b" << std::endl;

		++this->next_block_id;
		this->n_written_rcds += container.sizeWithoutCarryOver();

		l.unlock();
		cv_next_checkpoint.notify_all();
	}

	bool Write(containers::VariantBlock& block,
	           const containers::VcfContainer& container,
	           index::VariantIndexEntry& index_entry)
	{
		this->WriteBlock(block, index_entry); // write block
		this->UpdateIndex(container, index_entry); // Update index.
		return(this->writer->stream->good());
	}

	bool UpdateIndex(const containers::VcfContainer& container, index::VariantIndexEntry& index_entry){
		assert(this->writer != nullptr);
		index_entry.blockID         = this->writer->n_blocks_written;
		index_entry.byte_offset_end = this->writer->stream->tellp();
		index_entry.contigID        = container.front()->rid;
		index_entry.minPosition     = container.front()->pos;
		index_entry.n_variants      = container.sizeWithoutCarryOver();
		this->writer->index        += index_entry;

		index_entry.print(std::cerr);
		std::cerr << std::endl;

		index_entry.reset();
		++this->writer->n_blocks_written;
		this->writer->n_variants_written += container.sizeWithoutCarryOver();

		return true;
	}

	bool WriteBlock(containers::VariantBlock& block, index::VariantIndexEntry& index_entry){
		index_entry.byte_offset = this->writer->stream->tellp();
		block.write(*this->writer->stream);
		this->writer->WriteBlockFooter(block.footer_support);
		this->writer->WriteEndOfBlock(); // End-of-block marker.

		return(this->writer->stream->good());
	}

	bool WriteFinal(algorithm::VariantDigestManager& checksums){
		// Done importing
		this->writer->stream->flush();

		// Write global footer.
		core::Footer footer;
		footer.offset_end_of_data = this->writer->stream->tellp();
		footer.n_blocks           = this->writer->n_blocks_written;
		footer.n_variants         = this->writer->n_variants_written;
		assert(footer.n_blocks == this->writer->index.GetLinearSize());

		uint64_t last_pos = this->writer->stream->tellp();
		this->writer->WriteIndex(); // Write index.
		std::cerr << utility::timestamp("PROGRESS") << "Index size: " << utility::toPrettyDiskString((uint64_t)this->writer->stream->tellp() - last_pos) << "..." << std::endl;
		last_pos = this->writer->stream->tellp();
		checksums.finalize();       // Finalize SHA-512 digests.
		*this->writer->stream << checksums;
		std::cerr << utility::timestamp("PROGRESS") << "Checksum size: " << utility::toPrettyDiskString((uint64_t)this->writer->stream->tellp() - last_pos) << "..." << std::endl;
		last_pos = this->writer->stream->tellp();
		*this->writer->stream << footer; // Write global footer and EOF marker.
		std::cerr << utility::timestamp("PROGRESS") << "Footer size: " << utility::toPrettyDiskString((uint64_t)this->writer->stream->tellp() - last_pos) << "..." << std::endl;

		this->writer->stream->flush();
		return(this->writer->stream->good());
	}

public:
	uint64_t n_written_rcds;
	uint32_t next_block_id;
	std::atomic<bool> alive;
	std::mutex lock;
	std::condition_variable cv_next_checkpoint;

	VariantImportWriterInterface* writer; // writer

	// Stats
	support::VariantImporterContainerStats stats_basic;
	support::VariantImporterContainerStats stats_info;
	support::VariantImporterContainerStats stats_format;
};

struct yon_consumer_vcfc {
public:
	yon_consumer_vcfc(void) : n_rcds_processed(0), thread_id(0), b_shared(0), b_indiv(0), data_available(nullptr), data_pool(nullptr), global_header(nullptr), poolw(nullptr){}
	yon_consumer_vcfc(uint32_t thread_id, std::atomic<bool>& data_available, yon_pool_vcfc& data_pool) :
		n_rcds_processed(0), thread_id(thread_id), b_shared(0), b_indiv(0), data_available(&data_available), data_pool(&data_pool), global_header(nullptr),
		poolw(nullptr)
	{}

	yon_consumer_vcfc& operator+=(const yon_consumer_vcfc& other){
		this->n_rcds_processed += other.n_rcds_processed;
		this->b_shared += other.b_shared;
		this->b_indiv  += other.b_indiv;
		this->importer += other.importer;
		return(*this);
	}

	std::thread& Start(void){
		this->thread_ = std::thread(&yon_consumer_vcfc::Consume, this);
		return(this->thread_);
	}

private:
	bool Consume(void){
		algorithm::GenotypeSorter sorter;
		sorter.SetSamples(this->global_header->GetNumberSamples());

		while(data_available){
			yon_pool_vcfc_payload* d = data_pool->pop();
			if(d == nullptr){
				//std::cerr << "is exit conditon: nullptr" << std::endl;
				break;
			}
			assert(d->c != nullptr);
			uint64_t b_l_shared = 0, b_l_indiv = 0;
			for(int i = 0; i < d->c->sizeWithoutCarryOver(); ++i){
				b_l_shared += d->c->at(i)->shared.l;
				b_l_indiv  += d->c->at(i)->indiv.l;
			}
			this->b_shared += b_l_shared;
			this->b_indiv  += b_l_indiv;

			// Unpack bcf1_t records.
			n_rcds_processed += d->c->sizeWithoutCarryOver();
			for(int i = 0; i < d->c->sizeWithoutCarryOver(); ++i)
				bcf_unpack(d->c->at(i), BCF_UN_ALL);

			// Add data from the container to the local importer.
			assert(this->importer.Add(*d->c, d->block_id) == true);

			// Push
			this->poolw->cv_next_checkpoint.notify_all();
			this->poolw->emplace(d->block_id, *d->c, this->importer);

			// Cleanup data popped from the producer queue.
			delete d;
		}

		// Done.
		return true;
	}

public:
	uint64_t n_rcds_processed;
	uint32_t thread_id;
	uint64_t b_shared, b_indiv;
	std::atomic<bool>* data_available;
	yon_pool_vcfc* data_pool;
	std::thread thread_;
	io::VcfHeader* global_header;
	yon_writer_sync* poolw;
	VcfImporter importer;
};

} /* namespace Tachyon */

#endif /* VARIANT_IMPORTER_H_ */
