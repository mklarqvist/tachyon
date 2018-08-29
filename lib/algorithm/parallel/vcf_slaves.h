#ifndef ALGORITHM_PARALLEL_VCF_SLAVES_H_
#define ALGORITHM_PARALLEL_VCF_SLAVES_H_

#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "vcf_importer_slave.h"
#include "containers/vcf_container.h"
#include "core/variant_import_settings.h"
#include "core/variant_importer_container_stats.h"

namespace tachyon{

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
			// Exit condition to be triggered if the alive predicate
			// starts to evaluate as FALSE.
			if(this->n_c == 0 && this->alive == false) return true;
			return this->n_c != this->n_capacity;
		});

		// Exit condition when alive predicate evaluates as FALSE.
		if(this->n_c == 0 && this->alive == false){
			l.unlock();
			this->not_full.notify_all(); // flush
			this->not_empty.notify_all(); // flush
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
			// Exit condition to be triggered if the alive predicate
			// starts to evaluate as FALSE.
			if(this->n_c == 0 && this->alive == false) return true;
			return this->n_c != 0;
		});

		// Exit condition when alive predicate evaluates as FALSE.
		if(this->n_c == 0){
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

// Synchronised writer.
struct yon_writer_sync {
	yon_writer_sync() : n_written_rcds(0), next_block_id(0), alive(true), writer(nullptr){}
	~yon_writer_sync(){}

	void emplace(uint32_t block_id, const containers::VcfContainer& container, VcfImporterSlave& importer){
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

	VariantWriterInterface* writer; // writer

	// Stats
	support::VariantImporterContainerStats stats_basic;
	support::VariantImporterContainerStats stats_info;
	support::VariantImporterContainerStats stats_format;
};

/**<
 * Primary consumer instance for importing htslib bcf1_t records
 * into a tachyon archive. Invoking the Start() function will invoke
 * the spawning of a new thread that keeps popping yon_pool_vcfc_payload
 * objects from the shared pooled of payloads until end-of-file has been
 * reached in the production thread(s).
 */
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

			// Compute the number of bytes we have loaded in the
			// htslib bcf1_t records in the container.
			uint64_t b_l_shared = 0, b_l_indiv = 0;
			for(int i = 0; i < d->c->sizeWithoutCarryOver(); ++i){
				b_l_shared += d->c->at(i)->shared.l;
				b_l_indiv  += d->c->at(i)->indiv.l;
			}
			this->b_shared += b_l_shared;
			this->b_indiv  += b_l_indiv;

			// Unpack the htslib bcf1_t records.
			n_rcds_processed += d->c->sizeWithoutCarryOver();
			for(int i = 0; i < d->c->sizeWithoutCarryOver(); ++i)
				bcf_unpack(d->c->at(i), BCF_UN_ALL);

			// Add data from the container to the local importer.
			assert(this->importer.Add(*d->c, d->block_id) == true);
			this->poolw->cv_next_checkpoint.notify_all();

			// Push to the consumer pool. Notify all threads that are listening
			// that a new block has been placed with the given block_id. The thread
			// with the next block will push that to the writer.
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
	VcfImporterSlave importer;
};


}



#endif /* ALGORITHM_PARALLEL_VCF_SLAVES_H_ */
