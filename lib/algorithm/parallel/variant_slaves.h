#ifndef ALGORITHM_PARALLEL_VARIANT_SLAVES_H_
#define ALGORITHM_PARALLEL_VARIANT_SLAVES_H_

#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>

#include "variant_container.h"

namespace tachyon {

struct yon_pool_vblock_payload {
	yon_pool_vblock_payload() : block_id(0), c(nullptr) {}
	yon_pool_vblock_payload(const uint32_t bid, yon1_vb_t* vc) : block_id(bid), c(vc) {}
	~yon_pool_vblock_payload() { delete c; }
	yon_pool_vblock_payload(const yon_pool_vblock_payload& other) = delete; // copy is not allowed
	yon_pool_vblock_payload& operator=(const yon_pool_vblock_payload& other) = delete; // copy assign is not allowed
	yon_pool_vblock_payload(yon_pool_vblock_payload&& other) noexcept : block_id(other.block_id), c(nullptr) {
		std::swap(this->c, other.c);
	}
	yon_pool_vblock_payload& operator=(yon_pool_vblock_payload&& other) noexcept{
		if (this == &other) return(*this);

		this->block_id = other.block_id;
		delete this->c;
		this->c = nullptr;
		std::swap(this->c, other.c);
		return(*this);
	}

	uint32_t block_id;
	yon1_vb_t* c;
};

struct yon_pool_vblock {
public:
	yon_pool_vblock(void) :
		c(nullptr),
		n_capacity(0), n_c(0),
		front(0), rear(0), alive(false)
	{

	}

	yon_pool_vblock(uint32_t capacity) :
		c(new yon_pool_vblock_payload*[capacity]),
		n_capacity(capacity), n_c(0),
		front(0), rear(0), alive(false)
	{

	}

	yon_pool_vblock(const yon_pool_vblock& other) = delete; // disallow copy
	yon_pool_vblock& operator=(const yon_pool_vblock& other) = delete; // disallow assign copy

	yon_pool_vblock(yon_pool_vblock&& other) noexcept :
		c(nullptr),
		n_capacity(other.n_capacity), n_c(other.n_c),
		front(other.front), rear(other.rear), alive(other.alive.load())
	{
		std::swap(this->c, other.c);
	}

	yon_pool_vblock& operator=(yon_pool_vblock&& other) noexcept {
		if (this == &other) { return(*this); }
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

	~yon_pool_vblock() {
		delete [] this->c;
	}

	/**<
	 * Add payload to the cyclic queue. If the queue is full then
	 * wait until an item has been popped.
	 * @param data Input pointer to payload.
	 */
	void emplace(yon_pool_vblock_payload* data) {
		std::unique_lock<std::mutex> l(lock);

		this->not_full.wait(l, [this]() {
			// Exit condition to be triggered if the alive predicate
			// starts to evaluate as FALSE.
			if (this->n_c == 0 && this->alive == false) return true;
			return this->n_c != this->n_capacity;
		});

		// Exit condition when alive predicate evaluates as FALSE.
		if (this->n_c == 0 && this->alive == false) {
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
	 * Retrieve payload from the cyclic queue of shared resources.
	 * If the queue is empty * then wait until an item has been inserted.
	 * @return Returns a pointer to the retrieved payload or a nullpointer in the special case the producer operations has finished.
	 */
	yon_pool_vblock_payload* pop(void) {
		std::unique_lock<std::mutex> l(lock);

		this->not_empty.wait(l, [this]() {
			// Exit condition to be triggered if the alive predicate
			// starts to evaluate as FALSE.
			if (this->n_c == 0 && this->alive == false) return true;
			return this->n_c != 0;
		});

		// Exit condition when alive predicate evaluates as FALSE.
		if (this->n_c == 0) {
			l.unlock();
			this->not_empty.notify_all();
			this->not_full.notify_all();
			return nullptr;
		}

		//std::cerr << "popping" << std::endl;
		yon_pool_vblock_payload* result = this->c[this->front];
		this->c[this->front] = nullptr;
		this->front = (this->front + 1) % this->n_capacity;
		--this->n_c;

		l.unlock();
		this->not_full.notify_one();

		return result;
	}

public:
	yon_pool_vblock_payload** c;
	int n_capacity;
	int n_c;
	int front;
	int rear;
	std::atomic<bool> alive;
	std::mutex lock;
	std::condition_variable not_full;
	std::condition_variable not_empty;
};

struct yon_producer_vblock_interface {
public:
	yon_producer_vblock_interface(void) :
		all_finished(false),
		n_rcds_loaded(0),
		data_available(false),
		data_pool(0),
		dst_block_(nullptr)
	{}

	yon_producer_vblock_interface(uint32_t pool_size) :
		all_finished(false),
		n_rcds_loaded(0),
		data_available(false),
		data_pool(pool_size),
		dst_block_(nullptr)
	{}

	virtual ~yon_producer_vblock_interface() {}
	virtual std::thread& Start(void) =0;

public:
	bool all_finished;
	uint64_t n_rcds_loaded;
	std::atomic<bool> data_available;
	yon_pool_vblock data_pool;
	std::thread thread_;
	yon1_vb_t* dst_block_;
};

/**<
 * Single producer for producing raw (uncompressed and
 * encrypted) VcfContainer blocks. The producer push these
 * containers into the shared data pool of payloads.
 */
template <class T, class F = bool(T::*)(void)>
struct yon_producer_vblock : public yon_producer_vblock_interface {
public:
	typedef yon_producer_vblock self_type;
	typedef yon_producer_vblock_interface parent_type;

public:
	yon_producer_vblock(uint32_t pool_size) :
		parent_type(pool_size),
		instance_(nullptr)
	{}

	~yon_producer_vblock() {}

	std::thread& Start(void) {
		this->thread_ = std::thread(&yon_producer_vblock::Produce, this);
		return(this->thread_);
	}

	/**<
	 * Spawn worker reading data into the producer data pool.
	 * @return Returns a reference to the spawned thread.
	 */
	void Setup(F produce_function,
	           T& vreader,
	           yon1_vb_t& block)
	{
		this->instance_  = &vreader;
		this->func_      = produce_function;
		this->dst_block_ = &block;
	}

private:
	/**<
	 * Internal function that continues to produce raw vblockcontainer
	 * payloads to be inserted into the shared resource pool for
	 * consumers to retrieve. Continues until no more data is available.
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool Produce(void) {
		uint32_t n_blocks = 0;
		this->data_available  = true;
		this->data_pool.alive = true;
		while(true) {
			if ((*this->instance_.*this->func_)() == false) {
				// No more data is available or an error was seen.
				// Trigger shared resources flag alive to no longer
				// evaluate as true. This triggers the exit condition
				// in all consumers after the conditional variables
				// for empty and full are triggered.
				this->data_available  = false;
				this->data_pool.alive = false;
				// Keep flushing until the predicate for all_finished
				// evaluates TRUE. This occurs after all consumer threads
				// have been joined.
				while(data_pool.n_c && all_finished == false) {
					data_pool.not_empty.notify_all();
					data_pool.not_full.notify_all();
				}
				// All producer threads have been successfully joined.
				// It is now safe to exit this function and return to
				// the main thread that spawned this instance.
				break;
			}
			this->n_rcds_loaded += 0;
			// Peculiar syntax for adding a new payload to the cyclic queue. The
			// move semantics is required for intended functionality!
			//std::cerr << "emplacing: " << this->dst_block_->size() << " @ " <<  this->data_pool.n_c << std::endl;
			this->data_pool.emplace(new yon_pool_vblock_payload(n_blocks++, new yon1_vb_t(std::move(*this->dst_block_))));
			*this->dst_block_ = yon1_vb_t();
		}
		return true;
	}

public:
	F func_;
	T* instance_;
};

/**<
 * Primary consumer instance for importing htslib bcf1_t records
 * into a tachyon archive. Invoking the Start() function will invoke
 * the spawning of a new thread that keeps popping yon_pool_vcfc_payload
 * objects from the shared pooled of payloads until end-of-file has been
 * reached in the production thread(s).
 */
template <class T, class F = bool(T::*)(yon1_vb_t*&)>
struct yon_consumer_vblock {
public:
	yon_consumer_vblock(void) : n_rcds_processed(0), thread_id(0), data_available(nullptr), data_pool(nullptr), instance_(nullptr) {}
	yon_consumer_vblock(uint32_t thread_id, std::atomic<bool>& data_available, yon_pool_vblock& data_pool) :
		n_rcds_processed(0), thread_id(thread_id), data_available(&data_available), data_pool(&data_pool), instance_(nullptr)
	{}

	~yon_consumer_vblock() {}

	yon_consumer_vblock& operator+=(const yon_consumer_vblock& other) {
		this->n_rcds_processed += other.n_rcds_processed;
		return(*this);
	}

	std::thread& Start(F produce_function, T& vreader)
	{

		this->instance_  = &vreader;
		this->func_      = produce_function;

		this->thread_ = std::thread(&yon_consumer_vblock::Consume, this);
		return(this->thread_);
	}

private:
	/**<
	 * Internal function for consuming VcfContainers from the shared
	 * resource pool produced by the general producer(s). This function
	 * continually pops a payload from the shared pool of payloads and
	 * unpacks the htslib bcf1_t records and process that information
	 * into a valid Tachyon block. When this process is finished
	 * the resulting VariantBlock is pushed to the synchronized
	 * writer to be written to a byte stream.
	 * @return Returns TRUE upon success or FALSE oterhwise.
	 */
	bool Consume(void) {

		// Continue to consume data until there is no more
		// or an exit condition has been triggered.
		while(data_available) {
			// Pop a Vcf container payload from the shared
			// resource pool.
			yon_pool_vblock_payload* d = data_pool->pop();

			// If the payload is a nullpointer then an exit
			// condition has been triggered. If this is the
			// case we discontinue the consumption for this
			// consumer.
			if (d == nullptr) {
				break;
			}
			assert(d->c != nullptr);

			// Invoke consumer function of interest
			if ((*this->instance_.*this->func_)(d->c) == false) {
				std::cerr << "error occurred" << std::endl;
			}

			// Cleanup data popped from the producer queue.
			delete d;
		}

		// Done.
		return true;
	}

public:
	uint64_t n_rcds_processed;
	uint32_t thread_id;
	std::atomic<bool>* data_available;
	yon_pool_vblock* data_pool;
	std::thread thread_;
	F func_;
	T* instance_;
};

}



#endif /* ALGORITHM_PARALLEL_VARIANT_SLAVES_H_ */
