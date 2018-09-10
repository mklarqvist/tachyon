#ifndef CONTAINERS_VARIANT_CONTAINER_H_
#define CONTAINERS_VARIANT_CONTAINER_H_

#include "core/data_block_settings.h"
#include "core/variant_record.h"
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
	virtual int32_t GetSignedInt32NoConv(const uint32_t p) =0;
	virtual int32_t GetSignedInt32(const uint32_t p) =0;
	virtual int32_t GetInt32(const uint32_t p) =0;
	virtual float GetFloat(const uint32_t p) =0;

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
	inline int32_t GetSignedInt32NoConv(const uint32_t p){ return((int32_t)this->entries_[p]); }

	int32_t GetSignedInt32(const uint32_t p){
		if(sizeof(return_ptype) == sizeof(int32_t)){
			return(this->GetSignedInt32NoConv(p));
		}
		else {
			// If the data is missing in the native format.
			if((*this)[p] == std::numeric_limits<return_ptype>::min()){
				return(std::numeric_limits<int32_t>::min());
			}
			// If the data is EOV in the native format.
			else if((*this)[p] == std::numeric_limits<return_ptype>::min() + 1){
				return(std::numeric_limits<int32_t>::min() + 1);
			}
			else
				(*this)[p];
		}
	}

	inline int32_t GetInt32(const uint32_t p){ return(this->GetSignedInt32NoConv(p)); }
	inline float GetFloat(const uint32_t p){ return((*this)[p]); }

public:
	return_ptype* entries_;
};

class VariantContainer {
public:
	typedef containers::DataContainer dc_type;

public:
	VariantContainer() : n_variants_(0), n_capacity_(0), variants_(nullptr){}
	~VariantContainer(){
		delete [] variants_;
	}

	yon1_vnt_t& operator[](const uint32_t p){ return(this->variants_[p]); }
	const yon1_vnt_t& operator[](const uint32_t p) const{ return(this->variants_[p]); }

	void resize(const uint32_t new_size){
		if(new_size < this->n_capacity_){
			if(new_size < this->n_variants_){
				this->n_variants_ = new_size;
				return;
			}
			return;
		}

		yon1_vnt_t* temp = this->variants_;
		this->variants_ = new yon1_vnt_t[new_size];
		for(int i = 0; i < this->n_variants_; ++i)
			this->variants_[i] = std::move(temp[i]);

		delete [] temp;
		this->n_capacity_ = new_size;
	}

	bool Build(containers::VariantBlock& variant_block, const VariantHeader& header){
		// Interlace meta streams into the variant records.
		this->AddController(variant_block.base_containers[YON_BLK_CONTROLLER]);
		this->AddContigs(variant_block.base_containers[YON_BLK_CONTIG]);
		this->AddPositions(variant_block.base_containers[YON_BLK_POSITION]);
		this->AddQuality(variant_block.base_containers[YON_BLK_QUALITY]);
		this->AddFilterIds(variant_block.base_containers[YON_BLK_ID_FILTER]);
		this->AddFormatIds(variant_block.base_containers[YON_BLK_ID_FORMAT]);
		this->AddInfoIds(variant_block.base_containers[YON_BLK_ID_INFO]);
		this->AddPloidy(variant_block.base_containers[YON_BLK_GT_PLOIDY]);
		this->AddAlleles(variant_block.base_containers[YON_BLK_ALLELES]);
		this->AddNames(variant_block.base_containers[YON_BLK_NAMES]);
		this->AddRefAlt(variant_block.base_containers[YON_BLK_REFALT]);


		// Map local in block to local in variant.
		// Data is added to variants in the order they appear in the raw data block.
		// This is not necessarily corect as data may have been intended to be read in
		// a different order. We have to map from these global identifiers to the per-variant
		// local offsets by permuting the global to local order.
		std::vector< std::vector<int> > local_info_patterns(variant_block.load_settings->info_patterns_local.size());
		for(int i = 0; i < local_info_patterns.size(); ++i){
			std::vector<std::pair<int,int>> internal(variant_block.load_settings->info_patterns_local[i].size());
			for(int j = 0 ; j < internal.size(); ++j){
				internal[j] = std::pair<int,int>((*variant_block.footer.info_map)[variant_block.load_settings->info_patterns_local[i][j]], j);
			}
			std::sort(internal.begin(), internal.end());

			local_info_patterns[i].resize(internal.size());
			for(int j = 0; j < internal.size(); ++j){
				local_info_patterns[i][internal[j].second] = j;
			}
		}

		std::vector< std::vector<int> > local_format_patterns(variant_block.load_settings->format_patterns_local.size());
		for(int i = 0; i < local_format_patterns.size(); ++i){
			std::vector<std::pair<int,int>> internal(variant_block.load_settings->format_patterns_local[i].size());
			for(int j = 0 ; j < internal.size(); ++j){
				internal[j] = std::pair<int,int>((*variant_block.footer.format_map)[variant_block.load_settings->format_patterns_local[i][j]], j);
			}
			std::sort(internal.begin(), internal.end());

			local_format_patterns[i].resize(internal.size());
			for(int j = 0; j < internal.size(); ++j){
				local_format_patterns[i][internal[j].second] = j;
			}
		}

		// Allocate memory to support upcoming data to overload into the
		// containers. After memory allocation we provide the list of local
		// offsets that provide data.
		for(int i = 0; i < this->n_variants_; ++i){
			this->variants_[i].info   = new containers::PrimitiveContainerInterface*[variant_block.footer.n_info_streams];
			this->variants_[i].m_info = variant_block.footer.n_info_streams;
			this->variants_[i].fmt    = new containers::PrimitiveGroupContainerInterface*[variant_block.footer.n_format_streams];
			this->variants_[i].m_fmt  = variant_block.footer.n_format_streams;

			if(this->variants_[i].info_pid >= 0)
				this->variants_[i].info_ids = new std::vector<int>(local_info_patterns[this->variants_[i].info_pid]);

			if(this->variants_[i].fmt_pid >= 0)
				this->variants_[i].fmt_ids = new std::vector<int>(local_format_patterns[this->variants_[i].fmt_pid]);
		}

		this->AddFilter(variant_block, header);
		this->AddInfo(variant_block, header);
		this->AddFormat(variant_block, header);

		return true;
	}

	bool AddInfo(containers::VariantBlock& variant_block, const VariantHeader& header){
		for(int i = 0; i < variant_block.load_settings->info_id_local_loaded.size(); ++i){
			// Evaluate the set-membership of a given global key in the available Info patterns
			// described in the data container footer.
			std::vector<bool> matches = variant_block.InfoPatternSetMembership(variant_block.info_containers[i].header.GetGlobalKey());
			// Add Info data.
			this->AddInfoWrapper(variant_block.info_containers[variant_block.load_settings->info_id_local_loaded[i]], header, matches);
		}
		return true;
	}

	bool AddFilter(containers::VariantBlock& variant_block, const VariantHeader& header){
		for(int i = 0; i < this->n_variants_; ++i){
			const yon_blk_bv_pair& filter_patterns = variant_block.footer.filter_patterns[this->variants_[i].flt_pid];
			for(int j = 0; j < filter_patterns.pattern.size(); ++j){
				this->variants_[i].flt_hdr.push_back(&header.filter_fields_[filter_patterns.pattern[j]]);
			}
		}
		return true;
	}

	bool AddFormat(containers::VariantBlock& variant_block, const VariantHeader& header){
		for(int i = 0; i < variant_block.load_settings->format_id_local_loaded.size(); ++i){
			// If genotype data has been loaded.
			if(header.format_fields_[variant_block.format_containers[variant_block.load_settings->format_id_local_loaded[i]].header.GetGlobalKey()].id == "GT"){
				assert(variant_block.load_settings->loaded_genotypes);
				std::vector<bool> matches = variant_block.FormatPatternSetMembership(variant_block.format_containers[i].header.GetGlobalKey());

				// Add genotypes.
				this->AddGenotypes(variant_block, header);

				for(int j = 0; j < this->n_variants_; ++j){
					if(this->variants_[j].fmt_pid == -1){
						continue;
					} else if(matches[this->variants_[j].fmt_pid]){
						this->variants_[j].is_loaded_gt = true;
						this->variants_[j].fmt[this->variants_[j].n_fmt++] = new containers::PrimitiveGroupContainer<int32_t>();
						this->variants_[j].fmt_hdr.push_back(&header.format_fields_[variant_block.format_containers[variant_block.load_settings->format_id_local_loaded[i]].header.GetGlobalKey()]);

						// Check if genotype data has been loaded then checks
						// if this variant site has any genotypes set.
						if(variant_block.header.controller.has_gt_permuted)
							this->variants_[j].gt = this->variants_[i].gt_raw->GetObjects(*variant_block.gt_ppa);
						else
							this->variants_[j].gt = this->variants_[i].gt_raw->GetObjects(header.GetNumberSamples());

						this->variants_[j].gt->Evaluate();
						//records[i].gt->Expand();
					}
				}

				continue;
			}

			// Evaluate the set-membership of a given global key in the available Format patterns
			// described in the data container footer.
			std::vector<bool> matches = variant_block.FormatPatternSetMembership(variant_block.format_containers[i].header.GetGlobalKey());
			// Add Format data.
			this->AddFormatWrapper(variant_block.format_containers[variant_block.load_settings->format_id_local_loaded[i]], header, matches);
		}
		return true;
	}

	bool AddInfoWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
		if(container.data_uncompressed.size() == 0 &&
		   header.info_fields_[container.header.GetGlobalKey()].yon_type == YON_VCF_HEADER_FLAG)
		{
			for(int i = 0; i < this->n_variants_; ++i){
				if(this->variants_[i].info_pid == -1){
					continue;
				} else if(matches[this->variants_[i].info_pid]){
					this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<int32_t>(0);
					this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				}
			}
			return false;
		}

		if(container.data_uncompressed.size() == 0){
			return false;
		}

		if(container.header.data_header.HasMixedStride()){
			if(container.header.data_header.IsSigned()){
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->InfoSetup<int8_t>(container, header,  matches));  break;
				case(YON_TYPE_16B):    (this->InfoSetup<int16_t>(container, header,  matches));    break;
				case(YON_TYPE_32B):    (this->InfoSetup<int32_t>(container, header,  matches));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<int64_t>(container, header,  matches));    break;
				case(YON_TYPE_FLOAT):  (this->InfoSetup<float>(container, header,  matches));  break;
				case(YON_TYPE_DOUBLE): (this->InfoSetup<double>(container, header,  matches)); break;
				case(YON_TYPE_CHAR):   (this->InfoSetupString(container, header,  matches)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			} else {
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->InfoSetup<uint8_t>(container, header,  matches));   break;
				case(YON_TYPE_16B):    (this->InfoSetup<uint16_t>(container, header,  matches));    break;
				case(YON_TYPE_32B):    (this->InfoSetup<uint32_t>(container, header,  matches));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<uint64_t>(container, header,  matches));    break;
				case(YON_TYPE_FLOAT):  (this->InfoSetup<float>(container, header,  matches));  break;
				case(YON_TYPE_DOUBLE): (this->InfoSetup<double>(container, header,  matches)); break;
				case(YON_TYPE_CHAR):   (this->InfoSetupString(container, header,  matches)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			}
		} else {
			if(container.header.data_header.IsSigned()){
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->InfoSetup<int8_t>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_16B):    (this->InfoSetup<int16_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_32B):    (this->InfoSetup<int32_t>(container, header,  matches, container.header.data_header.stride));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<int64_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_FLOAT):  (this->InfoSetup<float>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_DOUBLE): (this->InfoSetup<double>(container, header,  matches, container.header.data_header.stride)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_CHAR):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			} else {
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->InfoSetup<uint8_t>(container, header,  matches, container.header.data_header.stride));   break;
				case(YON_TYPE_16B):    (this->InfoSetup<uint16_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_32B):    (this->InfoSetup<uint32_t>(container, header,  matches, container.header.data_header.stride));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<uint64_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_FLOAT):  (this->InfoSetup<float>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_DOUBLE): (this->InfoSetup<double>(container, header,  matches, container.header.data_header.stride)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_CHAR):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;

				}
			}
		}
		return true;
	}

	template <class return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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
		uint32_t i = 0;
		uint32_t stride_offset = 0;

		while(true){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, it->GetInt32(stride_offset));
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += it->GetInt32(stride_offset) * sizeof(return_ptype);
				++stride_offset;

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		delete it;
		return true;
	}

	template <class return_ptype>
	bool InfoSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size){
		if(container.strides_uncompressed.size() == 0)
			return false;

		uint32_t current_offset = 0;
		uint32_t i = 0;

		while(true){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<return_ptype>(container, current_offset, stride_size);
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size;

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		return true;
	}

	bool InfoSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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
		uint32_t i = 0;
		uint32_t stride_offset = 0;

		while(true){
			if(this->variants_[i].info_pid == -1){
				continue;
			} else if(matches[this->variants_[i].info_pid]){
				this->variants_[i].info[this->variants_[i].n_info++] = new containers::PrimitiveContainer<std::string>(&container.data_uncompressed[current_offset], it->GetInt32(stride_offset));
				this->variants_[i].info_hdr.push_back(&header.info_fields_[container.header.GetGlobalKey()]);
				current_offset += it->GetInt32(stride_offset) * sizeof(char);
				++stride_offset;

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		delete it;
		return true;
	}

	bool AddFormatWrapper(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
		if(container.data_uncompressed.size() == 0){
			return false;
		}

		if(container.header.data_header.HasMixedStride()){
			if(container.header.data_header.IsSigned()){
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->FormatSetup<int8_t>(container, header,  matches));  break;
				case(YON_TYPE_16B):    (this->FormatSetup<int16_t>(container, header,  matches));    break;
				case(YON_TYPE_32B):    (this->FormatSetup<int32_t>(container, header,  matches));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<int64_t>(container, header,  matches));    break;
				case(YON_TYPE_FLOAT):  (this->FormatSetup<float>(container, header,  matches));  break;
				case(YON_TYPE_DOUBLE): (this->FormatSetup<double>(container, header,  matches)); break;
				case(YON_TYPE_CHAR):   (this->FormatSetupString(container, header,  matches)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			} else {
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->FormatSetup<uint8_t>(container, header,  matches));   break;
				case(YON_TYPE_16B):    (this->FormatSetup<uint16_t>(container, header,  matches));    break;
				case(YON_TYPE_32B):    (this->FormatSetup<uint32_t>(container, header,  matches));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<uint64_t>(container, header,  matches));    break;
				case(YON_TYPE_FLOAT):  (this->FormatSetup<float>(container, header,  matches));  break;
				case(YON_TYPE_DOUBLE): (this->FormatSetup<double>(container, header,  matches)); break;
				case(YON_TYPE_CHAR):   (this->FormatSetupString(container, header,  matches)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			}
		} else {
			if(container.header.data_header.IsSigned()){
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->FormatSetup<int8_t>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_16B):    (this->FormatSetup<int16_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_32B):    (this->FormatSetup<int32_t>(container, header,  matches, container.header.data_header.stride));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<int64_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_FLOAT):  (this->FormatSetup<float>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_DOUBLE): (this->FormatSetup<double>(container, header,  matches, container.header.data_header.stride)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_CHAR):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;
				}
			} else {
				switch(container.header.data_header.GetPrimitiveType()){
				case(YON_TYPE_8B):     (this->FormatSetup<uint8_t>(container, header,  matches, container.header.data_header.stride));   break;
				case(YON_TYPE_16B):    (this->FormatSetup<uint16_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_32B):    (this->FormatSetup<uint32_t>(container, header,  matches, container.header.data_header.stride));    break;
				//case(YON_TYPE_64B):    (this->InfoSetup<uint64_t>(container, header,  matches, container.header.data_header.stride));    break;
				case(YON_TYPE_FLOAT):  (this->FormatSetup<float>(container, header,  matches, container.header.data_header.stride));  break;
				case(YON_TYPE_DOUBLE): (this->FormatSetup<double>(container, header,  matches, container.header.data_header.stride)); break;
				case(YON_TYPE_BOOLEAN):
				case(YON_TYPE_CHAR):
				case(YON_TYPE_STRUCT):
				case(YON_TYPE_UNKNOWN):
				default: std::cerr << "Disallowed type: " << (int)container.header.data_header.controller.type << std::endl; return false;

				}
			}
		}
		return true;
	}

	template <class return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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
		uint32_t i = 0;
		uint32_t stride_offset = 0;

		while(true){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), it->GetInt32(stride_offset));
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += it->GetInt32(stride_offset) * sizeof(return_ptype) * header.GetNumberSamples();
				++stride_offset;

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		delete it;
		return true;
	}

	template <class return_ptype>
	bool FormatSetup(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches, const uint32_t stride_size){
		if(container.strides_uncompressed.size() == 0)
			return false;

		uint32_t current_offset = 0;
		uint32_t i = 0;

		while(true){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<return_ptype>(container, current_offset, header.GetNumberSamples(), stride_size);
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += stride_size * header.GetNumberSamples();

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		return true;
	}

	bool AddGenotypes(containers::VariantBlock& block, const VariantHeader& header){
		const bool uniform_stride = block.base_containers[YON_BLK_GT_SUPPORT].header.data_header.IsUniform();
		containers::PrimitiveContainer<uint32_t> lengths(block.base_containers[YON_BLK_GT_SUPPORT]); // n_runs / objects size

		uint64_t offset_rle8     = 0; const char* const rle8     = block.base_containers[YON_BLK_GT_INT8].data_uncompressed.data();
		uint64_t offset_rle16    = 0; const char* const rle16    = block.base_containers[YON_BLK_GT_INT16].data_uncompressed.data();
		uint64_t offset_rle32    = 0; const char* const rle32    = block.base_containers[YON_BLK_GT_INT32].data_uncompressed.data();
		uint64_t offset_rle64    = 0; const char* const rle64    = block.base_containers[YON_BLK_GT_INT64].data_uncompressed.data();

		uint64_t offset_simple8  = 0; const char* const simple8  = block.base_containers[YON_BLK_GT_S_INT8].data_uncompressed.data();
		uint64_t offset_simple16 = 0; const char* const simple16 = block.base_containers[YON_BLK_GT_S_INT16].data_uncompressed.data();
		uint64_t offset_simple32 = 0; const char* const simple32 = block.base_containers[YON_BLK_GT_S_INT32].data_uncompressed.data();
		uint64_t offset_simple64 = 0; const char* const simple64 = block.base_containers[YON_BLK_GT_S_INT64].data_uncompressed.data();

		uint64_t offset_nploid8  = 0; const char* const nploid8  = block.base_containers[YON_BLK_GT_N_INT8].data_uncompressed.data();
		uint64_t offset_nploid16 = 0; const char* const nploid16 = block.base_containers[YON_BLK_GT_N_INT16].data_uncompressed.data();
		uint64_t offset_nploid32 = 0; const char* const nploid32 = block.base_containers[YON_BLK_GT_N_INT32].data_uncompressed.data();
		uint64_t offset_nploid64 = 0; const char* const nploid64 = block.base_containers[YON_BLK_GT_N_INT64].data_uncompressed.data();

		assert(block.base_containers[YON_BLK_GT_INT8].data_uncompressed.size()    % sizeof(uint8_t) == 0);
		assert(block.base_containers[YON_BLK_GT_INT16].data_uncompressed.size()   % sizeof(uint16_t)  == 0);
		assert(block.base_containers[YON_BLK_GT_INT32].data_uncompressed.size()   % sizeof(uint32_t)  == 0);
		assert(block.base_containers[YON_BLK_GT_INT64].data_uncompressed.size()   % sizeof(uint64_t)  == 0);
		assert(block.base_containers[YON_BLK_GT_S_INT8].data_uncompressed.size()  % sizeof(uint8_t) == 0);
		assert(block.base_containers[YON_BLK_GT_S_INT16].data_uncompressed.size() % sizeof(uint16_t)  == 0);
		assert(block.base_containers[YON_BLK_GT_S_INT32].data_uncompressed.size() % sizeof(uint32_t)  == 0);
		assert(block.base_containers[YON_BLK_GT_S_INT64].data_uncompressed.size() % sizeof(uint64_t)  == 0);

		uint64_t gt_offset = 0;
		uint8_t incrementor = 1;
		if(uniform_stride) incrementor = 0;

		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].controller.gt_available){
				// Case run-length encoding diploid and biallelic and no missing
				if(this->variants_[i].controller.gt_compression_type == TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_BIALLELIC){
					if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint8_t>( &rle8[offset_rle8], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_rle8 += lengths[gt_offset]*sizeof(uint8_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint16_t>( &rle16[offset_rle16], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_rle16 += lengths[gt_offset]*sizeof(uint16_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint32_t>( &rle32[offset_rle32], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_rle32 += lengths[gt_offset]*sizeof(uint32_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint64_t>( &rle64[offset_rle64], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_rle64 += lengths[gt_offset]*sizeof(uint64_t);
					} else {
						std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
						exit(1);
					}

				}
				// Case run-length encoding diploid and biallelic/EOV or n-allelic
				else if(this->variants_[i].controller.gt_compression_type == TACHYON_GT_ENCODING::YON_GT_RLE_DIPLOID_NALLELIC) {
					if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidSimple<uint8_t>( &simple8[offset_simple8], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple8 += lengths[gt_offset]*sizeof(uint8_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidSimple<uint16_t>( &simple16[offset_simple16], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
							offset_simple16 += lengths[gt_offset]*sizeof(uint16_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidSimple<uint32_t>( &simple32[offset_simple32], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple32 += lengths[gt_offset]*sizeof(uint32_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidSimple<uint64_t>( &simple64[offset_simple64], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple64 += lengths[gt_offset]*sizeof(uint64_t);
					} else {
						std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
						exit(1);
					}
				}
				// Case BCF-style encoding of diploids
				else if(this->variants_[i].controller.gt_compression_type == TACHYON_GT_ENCODING::YON_GT_BCF_DIPLOID) {
					if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidBCF<uint8_t>( &simple8[offset_simple8], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple8 += lengths[gt_offset]*sizeof(uint8_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidBCF<uint16_t>( &simple16[offset_simple16], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple16 += lengths[gt_offset]*sizeof(uint16_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidBCF<uint32_t>( &simple32[offset_simple32], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple32 += lengths[gt_offset]*sizeof(uint32_t);
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
						this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidBCF<uint64_t>( &simple64[offset_simple64], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_simple64 += lengths[gt_offset]*sizeof(uint64_t);
					}  else {
						std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
						exit(1);
					}
				}
				// Case RLE-encoding of nploids
				else if(this->variants_[i].controller.gt_compression_type == TACHYON_GT_ENCODING::YON_GT_RLE_NPLOID) {
					if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_BYTE){
						this->variants_[i].gt_raw = new containers::GenotypeContainerNploid<uint8_t>( &nploid8[offset_nploid8], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_nploid8 += lengths[gt_offset]*(sizeof(uint8_t) + this->variants_[i].n_base_ploidy*sizeof(uint8_t));
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U16){
						this->variants_[i].gt_raw = new containers::GenotypeContainerNploid<uint16_t>( &nploid16[offset_nploid16], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_nploid16 += lengths[gt_offset]*(sizeof(uint16_t) + this->variants_[i].n_base_ploidy*sizeof(uint8_t));
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U32){
						this->variants_[i].gt_raw = new containers::GenotypeContainerNploid<uint32_t>( &nploid32[offset_nploid32], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_nploid32 += lengths[gt_offset]*(sizeof(uint32_t) + this->variants_[i].n_base_ploidy*sizeof(uint8_t));
					} else if(this->variants_[i].controller.gt_primtive_type == TACHYON_GT_PRIMITIVE_TYPE::YON_GT_U64){
						this->variants_[i].gt_raw = new containers::GenotypeContainerNploid<uint64_t>( &nploid64[offset_nploid64], lengths[gt_offset], this->variants_[i].n_alleles, this->variants_[i].n_base_ploidy, this->variants_[i].controller.ToValue());
						offset_nploid64 += lengths[gt_offset]*(sizeof(uint64_t) + this->variants_[i].n_base_ploidy*sizeof(uint8_t));
					}  else {
						std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding primitive..." << std::endl;
						exit(1);
					}
				}
				// Case other potential encodings
				else {
					std::cerr << utility::timestamp("ERROR","GT") << "Unknown GT encoding family..." << std::endl;
					this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint8_t>( );
					exit(1);
				}

				// Increment offset
				gt_offset += incrementor;

			} else { // No GT available
				this->variants_[i].gt_raw = new containers::GenotypeContainerDiploidRLE<uint8_t>( );
			}
		}

		assert(offset_rle8     == block.base_containers[YON_BLK_GT_INT8].GetSizeUncompressed());
		assert(offset_rle16    == block.base_containers[YON_BLK_GT_INT16].GetSizeUncompressed());
		assert(offset_rle32    == block.base_containers[YON_BLK_GT_INT32].GetSizeUncompressed());
		assert(offset_rle64    == block.base_containers[YON_BLK_GT_INT64].GetSizeUncompressed());
		assert(offset_simple8  == block.base_containers[YON_BLK_GT_S_INT8].GetSizeUncompressed());
		assert(offset_simple16 == block.base_containers[YON_BLK_GT_S_INT16].GetSizeUncompressed());
		assert(offset_simple32 == block.base_containers[YON_BLK_GT_S_INT32].GetSizeUncompressed());
		assert(offset_simple64 == block.base_containers[YON_BLK_GT_S_INT64].GetSizeUncompressed());
		assert(offset_nploid8  == block.base_containers[YON_BLK_GT_N_INT8].GetSizeUncompressed());
		assert(offset_nploid16 == block.base_containers[YON_BLK_GT_N_INT16].GetSizeUncompressed());
		assert(offset_nploid32 == block.base_containers[YON_BLK_GT_N_INT32].GetSizeUncompressed());
		assert(offset_nploid64 == block.base_containers[YON_BLK_GT_N_INT64].GetSizeUncompressed());

		return true;
	}

	template <class T>
	bool AddBaseInteger(dc_type& container, void(yon1_vnt_t::*fnc)(const T v)){
		if(container.data_uncompressed.size() == 0)
			return false;

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

	bool FormatSetupString(dc_type& container, const VariantHeader& header, const std::vector<bool>& matches){
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
		uint32_t i = 0;
		uint32_t stride_offset = 0;

		while(true){
			if(this->variants_[i].fmt_pid == -1){
				continue;
			} else if(matches[this->variants_[i].fmt_pid]){
				this->variants_[i].fmt[this->variants_[i].n_fmt++] = new containers::PrimitiveGroupContainer<std::string>(container, current_offset, header.GetNumberSamples(), it->GetInt32(stride_offset));
				this->variants_[i].fmt_hdr.push_back(&header.format_fields_[container.header.GetGlobalKey()]);
				current_offset += it->GetInt32(stride_offset) * sizeof(char) * header.GetNumberSamples();
				++stride_offset;

				// Break condition
				if(current_offset == container.data_uncompressed.size())
					break;

				// Assertion of critical error
				assert(current_offset < container.data_uncompressed.size());
			}
			++i;
		}
		assert(current_offset == container.data_uncompressed.size());
		delete it;
		return true;
	}

	inline bool AddContigs(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetChromosome));
	}

	inline bool AddController(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetController));
	}

	inline bool AddPositions(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetPosition));
	}

	bool AddQuality(dc_type& container){
		if(container.data_uncompressed.size() == 0)
			return false;

		yon_cont_ref<float> it(container.data_uncompressed.data(), container.data_uncompressed.size());
		// If data is uniform.
		if(container.header.data_header.controller.uniform){
			it.is_uniform_ = true;

			for(int i = 0; i < this->n_variants_; ++i){
				this->variants_[i].qual = it[0];
			}
		} else {
			assert(this->n_variants_ == it.n_elements_);
			for(int i = 0; i < this->n_variants_; ++i){
				this->variants_[i].qual = it[i];
			}
		}

		return true;
	}

	bool AddRefAlt(dc_type& container){
		if(container.data_uncompressed.size() == 0)
			return false;

		yon_cont_ref<uint8_t> it(container.data_uncompressed.data(), container.data_uncompressed.size());

		uint32_t refalt_position = 0;
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].controller.alleles_packed){
				// load from special packed
				// this is always diploid
				this->variants_[i].n_alleles = 2;
				this->variants_[i].alleles   = new yon_allele[2];
				this->variants_[i].m_allele  = 2;

				// If data is <non_ref> or not
				if((it[refalt_position] & 15) != 5){
					//assert((refalt[refalt_position] & 15) < 5);
					const char ref = constants::REF_ALT_LOOKUP[it[refalt_position] & 15];
					this->variants_[i].alleles[1] = ref;

				} else {
					const std::string s = "<NON_REF>";
					this->variants_[i].alleles[1] = s;
				}

				// If data is <non_ref> or not
				if(((it[refalt_position] >> 4) & 15) != 5){
					//assert(((refalt[refalt_position] >> 4) & 15) < 5);
					const char alt = constants::REF_ALT_LOOKUP[(it[refalt_position] >> 4) & 15];
					this->variants_[i].alleles[0] = alt;
				} else {
					const std::string s = "<NON_REF>";
					this->variants_[i].alleles[1] = s;
				}
				// Do not increment in case this data is uniform
				if(it.is_uniform_ == false) ++refalt_position;
			}
			// otherwise load from literal cold
			else {
				// number of alleles is parsed from the stride container
			}
		}

		return true;
	}

	inline bool AddFilterIds(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetFilterPId));
	}

	inline bool AddFormatIds(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetFormatPId));
	}

	inline bool AddInfoIds(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetInfoPId));
	}

	inline bool AddPloidy(dc_type& container){
		return(this->AddBaseInteger(container, &yon1_vnt_t::SetBasePloidy));
	}

	bool AddAlleles(dc_type& container){
		if(container.data_uncompressed.size() == 0)
			return false;

		yon_cont_ref_iface* it = nullptr;
		switch(container.header.stride_header.controller.type){
		case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		}
		assert(it != nullptr);

		uint32_t offset = 0;
		uint32_t stride_offset = 0;
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			if(this->variants_[i].controller.alleles_packed == false){
				this->variants_[i].n_alleles = it->GetSignedInt32NoConv(stride_offset);
				this->variants_[i].alleles   = new yon_allele[it->GetSignedInt32NoConv(stride_offset)];
				this->variants_[i].m_allele  = it->GetSignedInt32NoConv(stride_offset);

				for(uint32_t j = 0; j < this->variants_[i].n_alleles; ++j){
					const uint16_t& l_string = *reinterpret_cast<const uint16_t* const>(&container.data_uncompressed[offset]);
					this->variants_[i].alleles[j].ParseFromBuffer(&container.data_uncompressed[offset]);
					offset += sizeof(uint16_t) + l_string;
				}
				++stride_offset;
			}
		}
		assert(offset == container.GetSizeUncompressed());
		delete it;
		return true;
	}

	bool AddNames(dc_type& container){
		if(container.data_uncompressed.size() == 0)
			return false;

		yon_cont_ref_iface* it = nullptr;
		switch(container.header.stride_header.controller.type){
		case(YON_TYPE_8B):  it = new yon_cont_ref<uint8_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_16B): it = new yon_cont_ref<uint16_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_32B): it = new yon_cont_ref<uint32_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		case(YON_TYPE_64B): it = new yon_cont_ref<uint64_t>(container.strides_uncompressed.data(), container.strides_uncompressed.size()); break;
		}
		assert(it != nullptr);

		uint32_t offset = 0;
		uint32_t stride_offset = 0;

		assert(it->n_elements_ == this->n_variants_);
		for(uint32_t i = 0; i < this->n_variants_; ++i){
			this->variants_[i].name = std::string(&container.data_uncompressed.data()[offset], it->GetSignedInt32NoConv(i));
			offset += it->GetSignedInt32NoConv(i);
		}
		assert(offset == container.data_uncompressed.size());
		delete it;
		return true;
	}

	/**< @brief Reads one or more separate digital objects from disk
	 * Primary function for reading partial data from disk. Data
	 * read in this way is not checked for integrity until later.
	 * @param stream   Input stream
	 * @param settings Settings record describing reading parameters
	 * @return         Returns FALSE if there was a problem, TRUE otherwise
	 */
	bool ReadBlock(std::ifstream& stream, const VariantHeader& header, DataBlockSettings& settings){
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

}

#endif /* CONTAINERS_VARIANT_CONTAINER_H_ */
