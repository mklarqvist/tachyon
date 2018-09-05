#ifndef ALGORITHM_PARALLEL_VARIANT_SLAVES_H_
#define ALGORITHM_PARALLEL_VARIANT_SLAVES_H_

#include <cstdint>
#include <thread>
#include <mutex>
#include <condition_variable>

#include "containers/variant_block.h"

namespace tachyon {

struct yon_pool_vblock_payload {
	yon_pool_vblock_payload() : block_id(0), c(nullptr){}
	yon_pool_vblock_payload(const uint32_t bid, containers::VariantBlock* vc) : block_id(bid), c(vc){}
	~yon_pool_vblock_payload(){ delete c; }
	yon_pool_vblock_payload(const yon_pool_vblock_payload& other) = delete; // copy is not allowed
	yon_pool_vblock_payload& operator=(const yon_pool_vblock_payload& other) = delete; // copy assign is not allowed
	yon_pool_vblock_payload(yon_pool_vblock_payload&& other) noexcept : block_id(other.block_id), c(nullptr){
		std::swap(this->c, other.c);
	}
	yon_pool_vblock_payload& operator=(yon_pool_vblock_payload&& other) noexcept{
		if(this == &other) return(*this);

		this->block_id = other.block_id;
		delete this->c;
		this->c = nullptr;
		std::swap(this->c, other.c);
		return(*this);
	}

	uint32_t block_id;
	containers::VariantBlock* c;
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

	~yon_pool_vblock(){
		delete [] this->c;
	}

	/**<
	 * Add payload to the cyclic queue. If the queue is full then
	 * wait until an item has been popped.
	 * @param data Input pointer to payload.
	 */
	void emplace(yon_pool_vblock_payload* data){
		std::cerr << "adding data: " << this->n_c << "/" << this->n_capacity << std::endl;
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
	 * Retrieve payload from the cyclic queue of shared resources.
	 * If the queue is empty * then wait until an item has been inserted.
	 * @return Returns a pointer to the retrieved payload or a nullpointer in the special case the producer operations has finished.
	 */
	yon_pool_vblock_payload* pop(void){
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

/**<
 * Single producer for producing raw (uncompressed and
 * encrypted) VcfContainer blocks. The producer push these
 * containers into the shared data pool of payloads.
 */
template <class T, class F = bool(T::*)(void)>
struct yon_producer_vblock {
public:
	typedef yon_producer_vblock self_type;

public:
	yon_producer_vblock(uint32_t pool_size) :
		all_finished(false),
		n_rcds_loaded(0),
		data_available(false),
		data_pool(pool_size),
		dst_block_(nullptr),
		instance_(nullptr)
	{}

	~yon_producer_vblock(){}

	/**<
	 * Spawn worker reading data into the producer data pool.
	 * @return Returns a reference to the spawned thread.
	 */
	std::thread& Start(F produce_function,
	                   T& vreader,
	                   containers::VariantBlock& block)
	{
		this->instance_  = &vreader;
		this->func_      = produce_function;
		this->dst_block_ = &block;

		this->thread_ = std::thread(&yon_producer_vblock::Produce, this);
		return(this->thread_);
	}

private:
	/**<
	 * Internal function that continues to produce raw vblockcontainer
	 * payloads to be inserted into the shared resource pool for
	 * consumers to retrieve. Continues until no more data is available.
	 * @return Returns TRUE if successful or FALSE otherwise.
	 */
	bool Produce(void){
		uint32_t n_blocks = 0;
		this->data_available  = true;
		this->data_pool.alive = true;
		while(true){
			// temp
			if(this->data_pool.n_c == this->data_pool.n_capacity){
				for(int i = 0; i < this->data_pool.n_c; ++i){
					delete this->data_pool.c[i];
					this->data_pool.c[i] = nullptr;
				}
				this->data_pool.n_c = 0;
			}

			if((*this->instance_.*this->func_)() == false){
				// temp
				for(int i = 0; i < this->data_pool.n_c; ++i){
					delete this->data_pool.c[i];
					this->data_pool.c[i] = nullptr;
				}
				this->data_pool.n_c = 0;
				std::cerr << "returning" << std::endl;
				// temp
				return false;

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
				while(data_pool.n_c && all_finished == false){
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
			std::cerr << "emplacing: " << this->dst_block_->size() << " @ " <<  this->data_pool.n_c << std::endl;
			this->data_pool.emplace(new yon_pool_vblock_payload(n_blocks++, new containers::VariantBlock(std::move(*this->dst_block_))));
			*this->dst_block_ = containers::VariantBlock();
		}
		return true;
	}

public:
	bool all_finished;
	uint64_t n_rcds_loaded;
	std::atomic<bool> data_available;
	yon_pool_vblock data_pool;
	std::thread thread_;
	F func_;
	T* instance_;
	containers::VariantBlock* dst_block_;
};

}



#endif /* ALGORITHM_PARALLEL_VARIANT_SLAVES_H_ */
