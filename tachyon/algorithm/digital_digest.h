#ifndef ALGORITHM_DigitalDigest_H_
#define ALGORITHM_DigitalDigest_H_

#include <cstring>
#include <openssl/sha.h>

#include "../io/basic_buffer.h"
#include "../containers/variantblock.h"

namespace tachyon{
namespace algorithm{

struct DigitalDigest{
private:
	typedef DigitalDigest             self_type;
	typedef containers::DataContainer container_type;
	typedef io::BasicBuffer           buffer_type;

public:
	DigitalDigest(void) : hasInitialized(true), hasFinished(false){ this->initialize(); }
	DigitalDigest(const self_type& other) :
		hasInitialized(other.hasInitialized),
		hasFinished(other.hasFinished),
		data_context(other.data_context),
		stride_context(other.stride_context)
	{
		memcpy(&this->data_digest[0],   other.data_digest,   64);
		memcpy(&this->stride_digest[0], other.stride_digest, 64);
	}

	DigitalDigest& operator=(const self_type& other){
		this->hasInitialized = other.hasInitialized;
		this->hasFinished = other.hasFinished;
		this->data_context = other.data_context;
		this->stride_context = other.stride_context;
		memcpy(&this->data_digest[0],   &other.data_digest[0],   64);
		memcpy(&this->stride_digest[0], &other.stride_digest[0], 64);
		return(*this);
	}

	~DigitalDigest(void){}

	/**<
	 * Initializes the SHA512 context. Calling this function
	 * is mandatory!
	 * @return
	 */
	inline bool initialize(){
		this->hasInitialized = true;

		if(!SHA512_Init(&this->data_context))   return false;
		if(!SHA512_Init(&this->stride_context)) return false;

		return true;
	}

	/**<
	 *
	 * @param data_buffer
	 * @param stride_buffer
	 * @param has_strides
	 * @return
	 */
	inline bool update(const buffer_type& data_buffer, const buffer_type& stride_buffer, const bool has_strides = true){
		if(!this->hasInitialized) this->initialize();

		if(!SHA512_Update(&this->data_context, (BYTE*)data_buffer.data(), data_buffer.size()))
			return false;

		if(has_strides){
			if(!SHA512_Update(&this->stride_context, (BYTE*)stride_buffer.data(), stride_buffer.size()))
				return false;
		}
		return true;
	}

	/**<
	 *
	 * @return
	 */
	inline bool finalize(){
		if(!this->hasInitialized) return false;
		if(this->hasFinished) return true;

		if(!SHA512_Final(&this->data_digest[0], &this->data_context))
			return false;

		if(!SHA512_Final(&this->stride_digest[0], &this->stride_context))
			return false;

		return true;
	}

	/**<
	 *
	 */
	inline void clear(void){
		this->hasFinished = false;
		this->finalize();
		memset(&this->data_digest[0],   0, 64);
		memset(&this->stride_digest[0], 0, 64);
		this->initialize();
		this->hasInitialized = true;
	}

private:
	/*
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << std::hex;
		for(U32 j = 0; j < 64; ++j)
			stream << (int)entry.data_digest[j];
		stream << std::dec << '\t' << std::hex;

		for(U32 j = 0; j < 64; ++j)
			stream << (int)entry.stride_digest[j];

		stream << std::hex;
		return(stream);
	}
	*/

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char* const>(&entry.data_digest),   64);
		stream.write(reinterpret_cast<const char* const>(&entry.stride_digest), 64);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.data_digest),   64);
		stream.read(reinterpret_cast<char*>(&entry.stride_digest), 64);
		return(stream);
	}

public:
	bool       hasInitialized;
	bool       hasFinished;
	SHA512_CTX data_context;
	SHA512_CTX stride_context;
	BYTE       data_digest[64];
	BYTE       stride_digest[64];
};

struct DigitalDigestPair{
private:
	typedef DigitalDigestPair         self_type;
	typedef DigitalDigest             digest_type;
	typedef containers::DataContainer container_type;

public:
	DigitalDigestPair(){}

	DigitalDigestPair(const self_type& other) :
		uncompressed(other.uncompressed),
		compressed(other.compressed)
	{

	}

	DigitalDigestPair& operator=(const self_type& other){
		this->uncompressed = other.uncompressed;
		this->compressed = other.compressed;
		return(*this);
	}

	~DigitalDigestPair(){}

	inline bool finalize(void){
		if(!this->uncompressed.finalize())
			return false;

		if(!this->compressed.finalize())
			return false;

		return true;
	}

	void operator+=(const container_type& container){
		this->compressed.update(container.buffer_data, container.buffer_strides, container.header.hasMixedStride());
		this->uncompressed.update(container.buffer_data_uncompressed, container.buffer_strides_uncompressed, container.header.hasMixedStride());
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream << entry.compressed;
		stream << entry.uncompressed;
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream >> entry.compressed;
		stream >> entry.uncompressed;
		return(stream);
	}


public:
	digest_type uncompressed;
	digest_type compressed;
};

class DigitalDigestManager{
private:
	typedef DigitalDigestManager self_type;

protected:
	typedef DigitalDigestPair value_type;
	typedef DigitalDigest     digest_type;
	typedef std::size_t       size_type;
	typedef value_type&       reference;
	typedef const value_type& const_reference;
	typedef value_type*       pointer;
	typedef const value_type* const_pointer;

public:
	DigitalDigestManager() :
		n_entries_(0),
		n_capacity_(100),
		__entries(new value_type[this->n_capacity_])
	{}

	DigitalDigestManager(const size_type start_capacity) :
		n_entries_(start_capacity),
		n_capacity_(start_capacity),
		__entries(new value_type[this->n_capacity_])
	{}

	DigitalDigestManager(const self_type& other) :
		n_entries_(other.n_entries_),
		n_capacity_(other.n_capacity_),
		__entries(new value_type[this->n_capacity_])
	{
		for(U32 i = 0; i < this->size(); ++i) this->__entries[i] = other.__entries[i];
	}

	virtual ~DigitalDigestManager(){ delete [] this->__entries; }

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	// Element access
	inline reference at(const size_type& position){ return(this->__entries[position]); }
	inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
	inline reference operator[](const size_type& position){ return(this->__entries[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
	inline pointer data(void){ return(this->__entries); }
	inline const_pointer data(void) const{ return(this->__entries); }
	inline reference front(void){ return(this->__entries[0]); }
	inline const_reference front(void) const{ return(this->__entries[0]); }
	inline reference back(void){ return(this->__entries[this->n_entries_ - 1]); }
	inline const_reference back(void) const{ return(this->__entries[this->n_entries_ - 1]); }

	// Capacity
	inline const bool empty(void) const{ return(this->n_entries_ == 0); }
	inline const size_type& size(void) const{ return(this->n_entries_); }
	inline const size_type& capacity(void) const{ return(this->n_capacity_); }

	// Iterator
	inline iterator begin(){ return iterator(&this->__entries[0]); }
	inline iterator end(){ return iterator(&this->__entries[this->n_entries_]); }
	inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries_]); }

	void finalize(void){
		for(U32 i = 0; i < this->size(); ++i) this->at(i).finalize();
	}

private:
	friend std::ostream& operator<<(std::ostream& out, const self_type& container){
		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_), sizeof(size_type));
		for(size_type i = 0; i < container.size(); ++i) out << container[i];

		return(out);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container){
		stream.read((char*)reinterpret_cast<size_type*>(&container.n_entries_), sizeof(size_type));
		delete [] container.__entries;
		container.__entries = new value_type[container.n_entries_];
		for(size_type i = 0; i < container.size(); ++i) stream >> container[i];

		return(stream);
	}

protected:
	size_type n_entries_;
	size_type n_capacity_;
	pointer   __entries;
};

class VariantDigitalDigestManager : public DigitalDigestManager{
private:
	typedef VariantDigitalDigestManager self_type;
	typedef DigitalDigestManager        parent_type;
	typedef containers::VariantBlock    variant_block_type;

public:
	VariantDigitalDigestManager() :
		parent_type(100),
		n_entries_info_(0),
		n_entries_format_(0),
		n_capacity_info_(100),
		n_capacity_format(100),
		__entries_info(new value_type[this->n_entries_info_]),
		__entries_format(new value_type[this->n_entries_format_])
	{
	}

	VariantDigitalDigestManager(const size_type base_capacity) :
		parent_type(base_capacity),
		n_entries_info_(base_capacity),
		n_entries_format_(base_capacity),
		n_capacity_info_(base_capacity),
		n_capacity_format(base_capacity),
		__entries_info(new value_type[this->n_entries_info_]),
		__entries_format(new value_type[this->n_entries_format_])
	{
	}

	VariantDigitalDigestManager(const size_type base_capacity, const size_type capacity_info, const size_type capacity_format) :
		parent_type(base_capacity),
		n_entries_info_(capacity_info),
		n_entries_format_(capacity_format),
		n_capacity_info_(capacity_info),
		n_capacity_format(capacity_format),
		__entries_info(new value_type[this->n_entries_info_]),
		__entries_format(new value_type[this->n_entries_format_])
	{
	}

	VariantDigitalDigestManager(const self_type& other) :
		parent_type(other),
		n_entries_info_(other.n_entries_info_),
		n_entries_format_(other.n_entries_format_),
		n_capacity_info_(other.n_capacity_info_),
		n_capacity_format(other.n_capacity_format),
		__entries_info(new value_type[this->n_entries_info_]),
		__entries_format(new value_type[this->n_entries_format_])
	{
		for(U32 i = 0; i < this->n_entries_info_; ++i) this->__entries_info[i] = other.__entries_info[i];
		for(U32 i = 0; i < this->n_entries_format_; ++i) this->__entries_format[i] = other.__entries_format[i];
	}

	~VariantDigitalDigestManager(){
		delete [] this->__entries_info;
		delete [] this->__entries_format;
	}

	void finalize(void){
		parent_type::finalize();
		for(U32 i = 0; i < this->n_capacity_info_; ++i) this->atINFO(i).finalize();
		for(U32 i = 0; i < this->n_capacity_format ; ++i) this->atFORMAT(i).finalize();
	}

	inline const_reference atINFO(const U32 position) const{ return(this->__entries_info[position]); }
	inline const_reference atFORMAT(const U32 position) const{ return(this->__entries_format[position]); }
	inline reference atINFO(const U32 position){ return(this->__entries_info[position]); }
	inline reference atFORMAT(const U32 position){ return(this->__entries_format[position]); }

	void operator+=(const variant_block_type& block){
		this->at(1)  += block.meta_contig_container;
		this->at(2)  += block.meta_positions_container;
		this->at(3)  += block.meta_names_container;
		this->at(4)  += block.meta_refalt_container;
		this->at(5)  += block.meta_controller_container;
		this->at(6)  += block.meta_quality_container;
		this->at(7)  += block.meta_names_container;
		this->at(8)  += block.meta_alleles_container;
		this->at(9)  += block.meta_info_map_ids;
		this->at(10) += block.meta_format_map_ids;
		this->at(11) += block.meta_filter_map_ids;
		this->at(12) += block.gt_support_data_container;
		this->at(13) += block.gt_rle8_container;
		this->at(14) += block.gt_rle16_container;
		this->at(15) += block.gt_rle32_container;
		this->at(16) += block.gt_rle64_container;
		this->at(17) += block.gt_simple8_container;
		this->at(18) += block.gt_simple16_container;
		this->at(19) += block.gt_simple32_container;
		this->at(20) += block.gt_simple64_container;

		for(U32 i = 0; i < block.footer.n_info_streams; ++i) this->__entries_info[block.footer.info_offsets[i].data_header.global_key] += block.info_containers[i];
		for(U32 i = 0; i < block.footer.n_format_streams; ++i) this->__entries_format[block.footer.format_offsets[i].data_header.global_key] += block.format_containers[i];
	}

	friend std::ostream& operator<<(std::ostream& out, const self_type& container){
		const parent_type* const parent = reinterpret_cast<const parent_type* const>(&container);
		out << *parent;

		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_info_), sizeof(size_type));
		out.write((const char* const)reinterpret_cast<const size_type* const>(&container.n_entries_format_), sizeof(size_type));
		for(size_type i = 0; i < container.n_entries_info_; ++i)   out << container.atINFO(i);
		for(size_type i = 0; i < container.n_entries_format_; ++i) out << container.atFORMAT(i);

		return(out);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& container){
		parent_type* parent = reinterpret_cast<parent_type*>(&container);
		stream >> *parent;

		stream.read((char*)reinterpret_cast<size_type*>(&container.n_entries_info_), sizeof(size_type));
		stream.read((char*)reinterpret_cast<size_type*>(&container.n_capacity_format), sizeof(size_type));

		delete [] container.__entries_info;
		delete [] container.__entries_format;
		container.__entries_info = new value_type[container.n_entries_info_];
		container.__entries_format = new value_type[container.n_entries_format_];
		for(size_type i = 0; i < container.n_entries_info_; ++i)   stream >> container.atINFO(i);
		for(size_type i = 0; i < container.n_entries_format_; ++i) stream >> container.atFORMAT(i);

		return(stream);
	}

private:
	size_type n_entries_info_;
	size_type n_entries_format_;
	size_type n_capacity_info_;
	size_type n_capacity_format;
	pointer   __entries_info;
	pointer   __entries_format;
};


}
}

#endif /* ALGORITHM_DigitalDigest_H_ */
