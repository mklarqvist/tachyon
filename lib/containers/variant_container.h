#ifndef CONTAINERS_VARIANT_CONTAINER_H_
#define CONTAINERS_VARIANT_CONTAINER_H_

#include "core/data_block_settings.h"
#include "core/variant_record.h"
#include "containers/variant_block.h"
#include "containers/genotype_container_diploid_bcf.h"
#include "containers/genotype_container_diploid_simple.h"
#include "containers/genotype_container_diploid_rle.h"
#include "containers/genotype_container_nploid.h"

#include "primitive_group_container_string.h"

namespace tachyon{

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
	yon_cont_ref() : entries_(nullptr){}
	yon_cont_ref(char* data, uint64_t l_data) :
		yon_cont_ref_iface(data, l_data),
		entries_(reinterpret_cast<return_ptype*>(this->data_))
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

	inline return_ptype& operator[](const uint32_t p){ return(this->entries_[p]); }
	inline const return_ptype& operator[](const uint32_t p) const{ return(this->entries_[p]); }
	inline int32_t GetInt32(const uint32_t p){ return((int32_t)this->entries_[p]); }

public:
	return_ptype* entries_;
};

class VariantContainer {
public:
	typedef containers::DataContainer dc_type;

public:
	VariantContainer() : n_variants_(0), n_capacity_(0), variants_(nullptr){}
	VariantContainer(const uint32_t size) : n_variants_(size), n_capacity_(size), variants_(new yon1_vnt_t[size]){}
	~VariantContainer(){
		delete [] variants_;
	}

	const size_t size(void) const { return(this->n_variants_); }
	const size_t capacity(void) const { return(this->n_capacity_); }

	yon1_vnt_t& operator[](const uint32_t p){ return(this->variants_[p]); }
	const yon1_vnt_t& operator[](const uint32_t p) const{ return(this->variants_[p]); }

	void resize(const uint32_t new_size);

	bool Build(containers::VariantBlock& variant_block, const VariantHeader& header);
	bool PermuteOrder(const containers::VariantBlock& variant_block);

	bool AddInfo(containers::VariantBlock& variant_block, const VariantHeader& header);
	bool AddFilter(containers::VariantBlock& variant_block, const VariantHeader& header);
	bool AddFormat(containers::VariantBlock& variant_block, const VariantHeader& header);
	bool AddInfoWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool AddFormatWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool InfoSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool InfoSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);

	template <class return_ptype, class intrinsic_ptype = return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size);

	bool AddGenotypes(containers::VariantBlock& block, const VariantHeader& header);

	template <class T>
	bool AddBaseInteger(dc_type& container, void(yon1_vnt_t::*fnc)(const T v)){
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

			for(int i = 0; i < this->n_variants_; ++i){
				(variants_[i].*fnc)(it->GetInt32(0));
			}
		} else {
			assert(this->n_variants_ == it->n_elements_);
			for(int i = 0; i < this->n_variants_; ++i){
				(variants_[i].*fnc)(it->GetInt32(i));
			}
		}
		delete it;
		return true;
	}

	bool FormatSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches);
	bool FormatSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride);

	inline bool AddContigs(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetChromosome)); }
	inline bool AddController(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetController)); }
	inline bool AddPositions(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetPosition)); }

	bool AddQuality(dc_type& container);
	bool AddRefAlt(dc_type& container);

	inline bool AddFilterIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFilterPId)); }
	inline bool AddFormatIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetFormatPId)); }
	inline bool AddInfoIds(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetInfoPId)); }
	inline bool AddPloidy(dc_type& container){ return(this->AddBaseInteger(container, &yon1_vnt_t::SetBasePloidy)); }

	bool AddAlleles(dc_type& container);
	bool AddNames(dc_type& container);

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity until later.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	inline bool ReadBlock(std::ifstream& stream, const VariantHeader& header, DataBlockSettings& settings){
		return(this->block_.read(stream, settings, header));
	}


public:
	containers::VariantBlock block_;
	// External memory allocation for linear use of lazy-evaluated
	// expansion of genotype records. This is critical when the sample
	// numbers are becoming large as allocating/deallocating hundreds
	// of thousands of pointers for every variant is very time consuming.
	//yon_gt_rcd** gt_exp;

	size_t n_variants_;
	size_t n_capacity_;
	yon1_vnt_t* variants_;
};

template <class return_ptype, class intrinsic_ptype>
bool VariantContainer::InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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

	for(int i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].info_pid == -1){
			continue;
		} else if(matches[this->variants_[i].info_pid]){
			this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, it->GetInt32(stride_offset));
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
bool VariantContainer::InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size){
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, 0, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype);
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

template <class return_ptype, class intrinsic_ptype>
bool VariantContainer::FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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

	for(int i = 0; i < this->n_variants_; ++i){
		if(this->variants_[i].fmt_pid == -1){
			continue;
		} else if(matches[this->variants_[i].fmt_pid]){
			this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), it->GetInt32(stride_offset));
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
bool VariantContainer::FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size){
	uint32_t current_offset = 0;

	if(container.header.data_header.IsUniform()){
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, 0, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
			}
		}
	} else {
		for(int i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * sizeof(intrinsic_ptype) * header.GetNumberSamples();
			}
		}
		assert(current_offset == container.data_uncompressed.size());
	}
	return true;
}

}

#endif /* CONTAINERS_VARIANT_CONTAINER_H_ */
