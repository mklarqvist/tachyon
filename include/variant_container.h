/*
Copyright (C) 2017-current Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#ifndef TACHYON_VARIANT_CONTAINER_H_
#define TACHYON_VARIANT_CONTAINER_H_

#include <vector>
#include <string>
#include <cstring>
#include <cassert>
#include <unordered_map>
#include <cmath>
#include <fstream>

#include "data_container.h"
#include "variant_block.h"
#include "primitive_container.h"
#include "variant_record.h"

namespace tachyon {

//
class yon_cont_ref_iface {
public:
	yon_cont_ref_iface() : is_uniform_(false), n_offset_(0), n_elements_(0), b_data_(0), data_(nullptr){}
	yon_cont_ref_iface(char* data, uint64_t l_data) : is_uniform_(false), n_offset_(0), n_elements_(0), b_data_(l_data), data_(data){}
	virtual ~yon_cont_ref_iface(){
		// do not delete data. it is not owned by this
	}

	// next()
	// -> stride of response + response pointer
	virtual bool next() =0;
	virtual int32_t GetInt32(const uint32_t p) =0;

public:
	bool is_uniform_;
	uint32_t n_offset_;
	uint32_t n_elements_;
	uint64_t b_data_;
	char* data_;
};

template <class return_ptype>
class yon_cont_ref : public yon_cont_ref_iface {
public:
	yon_cont_ref() : variants_(nullptr){}
	yon_cont_ref(char* data, uint64_t l_data) :
		yon_cont_ref_iface(data, l_data),
		variants_(reinterpret_cast<return_ptype*>(this->data_))
	{
		this->n_elements_ = l_data / sizeof(return_ptype);
		assert(l_data % sizeof(return_ptype) == 0);
	}
	~yon_cont_ref(){}

	bool next(){
		if(this->is_uniform_) return true;
		if(this->n_offset_ == this->n_elements_) return false;
		++this->n_offset_;
		return true;
	}

	inline return_ptype& operator[](const uint32_t p){ return(this->variants_[p]); }
	inline const return_ptype& operator[](const uint32_t p) const{ return(this->variants_[p]); }
	inline int32_t GetInt32(const uint32_t p){ return((int32_t)this->variants_[p]); }

public:
	return_ptype* variants_;
};

class yon1_vc_t {
public:
	typedef yon1_dc_t         dc_type;
	typedef yon1_vc_t         self_type;
    typedef std::size_t       size_type;
    typedef yon1_vnt_t        value_type;
    typedef value_type&       reference;
    typedef const value_type& const_reference;
    typedef value_type*       pointer;
    typedef const value_type* const_pointer;
    typedef std::ptrdiff_t    difference_type;

    typedef yonRawIterator<value_type>       iterator;
	typedef yonRawIterator<const value_type> const_iterator;

public:
	yon1_vc_t() : n_variants_(0), n_capacity_(0), variants_(nullptr){}
	yon1_vc_t(const uint32_t size) :
		n_variants_(0), n_capacity_(size),
		//variants_(new yon1_vnt_t[size])
		variants_(static_cast<yon1_vnt_t*>(::operator new[](this->n_capacity_*sizeof(yon1_vnt_t))))
	{

	}

	yon1_vc_t(yon1_vb_t& variant_block, const yon_vnt_hdr_t& header) :
		n_variants_(variant_block.header.n_variants), n_capacity_(variant_block.header.n_variants + 64),
		//variants_(new yon1_vnt_t[size])
		variants_(static_cast<yon1_vnt_t*>(::operator new[](this->n_capacity_*sizeof(yon1_vnt_t))))
	{
		if(this->n_variants_){
			for(uint32_t i = 0; i < this->n_variants_; ++i)
				new( &this->variants_[i] ) yon1_vnt_t( );

			this->Build(variant_block, header);
		}
	}

	~yon1_vc_t(){
		//delete [] variants_;
		for(std::size_t i = 0; i < this->size(); ++i)
			((this->variants_ + i)->~yon1_vnt_t)();

		::operator delete[](static_cast<void*>(this->variants_));
	}

	yon1_vc_t& operator+=(const yon1_vnt_t& rec){
		if(this->variants_ == nullptr || this->capacity() == 0)
			this->reserve();

		if(this->size() == this->capacity())
			this->reserve();

		new( &this->variants_[this->n_variants_] ) yon1_vnt_t( rec );
		++this->n_variants_;

		return(*this);
	}

	const size_t size(void) const { return(this->n_variants_); }
	const size_t capacity(void) const { return(this->n_capacity_); }

	void clear(void){
		for(std::size_t i = 0; i < this->size(); ++i)
			((this->variants_ + i)->~yon1_vnt_t)();

		this->n_variants_ = 0;
		this->block_.clear();
	}

	inline void reserve(void){ this->reserve(n_capacity_ + 1000); }
	void reserve(const uint32_t new_size);
	void resize(const uint32_t new_size);

	bool Build(yon1_vb_t& variant_block, const yon_vnt_hdr_t& header);

	/**< @brief Reads a target Tachyon block from disk.
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity in this function.
	 * @param stream    Src stream handle.
	 * @param header    Src global header.
	 * @param settings  Settings record describing reading parameters.
	 * @return          Returns TRUE upon success or FALSE otherwise.
	 */
	inline bool ReadBlock(std::ifstream& stream, const yon_vnt_hdr_t& header, yon_vb_settings& settings){
		return(this->block_.read(stream, settings, header));
	}

	/**<
	 * Prepares the data block structure in the variant container for writing.
	 * Internally permutes genotypes, perform compression, encryption, and prepares
	 * the internal header and footer structures.
	 * @param n_samples         Number of samples in the file.
	 * @param compression_level Dst compression level of data block.
	 * @return                  Returns TRUE upon success or FALSE otherwise.
	 */
	bool PrepareWritableBlock(const uint32_t n_samples, const uint32_t compression_level);

	 // Element access
	inline reference at(const size_type& position){ return(this->variants_[position]); }
	inline const_reference at(const size_type& position) const{ return(this->variants_[position]); }
	inline reference operator[](const size_type& position){ return(this->variants_[position]); }
	inline const_reference operator[](const size_type& position) const{ return(this->variants_[position]); }
	inline pointer data(void){ return(this->variants_); }
	inline const_pointer data(void) const{ return(this->variants_); }
	inline reference front(void){ return(this->variants_[0]); }
	inline const_reference front(void) const{ return(this->variants_[0]); }
	inline reference back(void){ return(this->variants_[(n_variants_ == 0 ? 0 : n_variants_ - 1)]); }
	inline const_reference back(void) const{ return(this->variants_[(n_variants_ == 0 ? 0 : n_variants_ - 1)]); }

	// Iterator
	inline iterator begin(){ return iterator(&this->variants_[0]); }
	inline iterator end(){ return iterator(&this->variants_[this->n_variants_]); }
	inline const_iterator begin() const{ return const_iterator(&this->variants_[0]); }
	inline const_iterator end() const{ return const_iterator(&this->variants_[this->n_variants_]); }
	inline const_iterator cbegin() const{ return const_iterator(&this->variants_[0]); }
	inline const_iterator cend() const{ return const_iterator(&this->variants_[this->n_variants_]); }

private:
	bool PermuteOrder(const yon1_vb_t& variant_block);

	bool AddInfo(yon1_vb_t& variant_block, const yon_vnt_hdr_t& header);
	bool AddFilter(yon1_vb_t& variant_block, const yon_vnt_hdr_t& header);
	bool AddFormat(yon1_vb_t& variant_block, const yon_vnt_hdr_t& header);
	bool AddInfoWrapper(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);
	bool AddFormatWrapper(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool InfoSetupString(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);
	bool InfoSetupString(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches, const uint32_t stride);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool AddGenotypes(yon1_vb_t& block, const yon_vnt_hdr_t& header);

	template <class T>
	bool AddBaseInteger(dc_type& container, void(yon1_vnt_t::*fnc)(const T v));

	bool FormatSetupString(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches);
	bool FormatSetupString(dc_type& container, const yon_vnt_hdr_t& header, const std::vector<bool>& matches, const uint32_t stride);

	inline bool AddContigs(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetChromosome)); }
	inline bool AddController(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetController)); }
	inline bool AddPositions(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetPosition)); }

	bool AddQuality(dc_type& container);
	bool AddRefAlt(dc_type& container);

	inline bool AddFilterIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFilterPatternId)); }
	inline bool AddFormatIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFormatPatternId)); }
	inline bool AddInfoIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetInfoPatternId)); }
	inline bool AddPloidy(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetBasePloidy)); }

	bool AddAlleles(dc_type& container);
	bool AddNames(dc_type& container);

public:
	yon1_vb_t block_;
	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	//yon_gt_rcd** gt_exp;

	size_t n_variants_;
	size_t n_capacity_;
	yon1_vnt_t* variants_;
};


// IMPLEMENTATION -------------------------------------------------------------


template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::InfoSetup(dc_type& container,
                                 const yon_vnt_hdr_t& header,
                                 const std::vector<bool>& matches)
{
	if(container.strides_uncompressed.size() == 0)
		return false;

	yon_cont_ref_iface* it = nullptr;
	switch(container.header.stride_header.controller.type){
	case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	}
	assert(it != nullptr);

	uint32_t current_offset = 0;
	uint32_t stride_offset = 0;

	for(uint32_t i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].info_pid == -1){
			continue;
		} else if(matches[this->variants_[i].info_pid]){
			this->variants_[i].info[this->variants_[i].n_info++] = new PrimitiveContainer<return_ptype>(container, current_offset, it->GetInt32(stride_offset));
			this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
			current_offset += it->GetInt32(stride_offset) * sizeof(intrinsic_ptype);
			++stride_offset;
		}
	}
	assert(current_offset == container.data_uncompressed.size());
	delete it;
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::InfoSetup(dc_type& container,
                                 const yon_vnt_hdr_t& header,
                                 const std::vector<bool>& matches,
                                 const uint32_t stride_size)
{
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new PrimitiveContainer<return_ptype>(container, 0, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new PrimitiveContainer<return_ptype>(container, current_offset, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype);
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::FormatSetup(dc_type& container,
                                   const yon_vnt_hdr_t& header,
                                   const std::vector<bool>& matches)
{
	if(container.strides_uncompressed.size() == 0)
		return false;

	yon_cont_ref_iface* it = nullptr;
	switch(container.header.stride_header.controller.type){
	case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
	}
	assert(it != nullptr);

	uint32_t current_offset = 0;
	uint32_t stride_offset = 0;

	for(uint32_t i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].fmt_pid == -1){
			continue;
		} else if(matches[this->variants_[i].fmt_pid]){
			this->variants_[i].fmt[this->variants_[i].n_fmt++] = new PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), it->GetInt32(stride_offset));
			this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
			current_offset += it->GetInt32(stride_offset) * sizeof(intrinsic_ptype) * header.GetNumberSamples();
			++stride_offset;
		}
	}
	assert(current_offset == container.data_uncompressed.size());
	delete it;
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool yon1_vc_t::FormatSetup(dc_type& container,
                            const yon_vnt_hdr_t& header,
                            const std::vector<bool>& matches,
                            const uint32_t stride_size)
{
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new PrimitiveGroupContainer<return_ptype>(container, 0, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype) * header.GetNumberSamples();
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

template <class T>
bool yon1_vc_t::AddBaseInteger(dc_type& container, void(yon1_vnt_t::*fnc)(const T v)){
	if(container.data_uncompressed.size() == 0)
		return false;

	assert(container.header.HasMixedStride() == false);

	yon_cont_ref_iface* it = nullptr;
	switch(container.header.data_header.controller.type){
	case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
	case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
	case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
	case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.data_uncompressed.data(), container.data_uncompressed.size()); break;
	}
	assert(it != nullptr);

	// If data is uniform.
	if(container.header.data_header.controller.uniform){
		it->is_uniform_ = true;

		for(uint32_t i = 0; i < this->n_variants_; ++i){
			(variants_[i].*fnc)(it->GetInt32(0));
		}
	} else {
		assert(this->n_variants_ == it->n_elements_);
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			(variants_[i].*fnc)(it->GetInt32(i));
		}
	}
	delete it;
	return true;
}

}



#endif /* CONTAINERS_VARIANT_CONTAINER_H_ */
