#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/timer.h"
#include "containers/abstract_integer_container.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/meta_cold_container.h"
#include "containers/meta_hot_container.h"
#include "containers/primitive_group_container.h"
#include "containers/meta_container.h"
#include "containers/checksum_container.h"
#include "containers/variantblock.h"
#include "core/genotype_object.h"
#include "core/header/tachyon_header.h"
#include "core/footer/footer.h"
#include "math/fisher.h"
#include "math/square_matrix.h"
#include "math/basic_vector_math.h"
#include "utility/support_vcf.h"
#include "iterator/IteratorIntegerReference.h"
#include "index/index.h"

namespace tachyon{

class TachyonReader{
	typedef TachyonReader                       self_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef io::BasicBuffer                     buffer_type;
	typedef core::TachyonHeader                 header_type;
	typedef core::Footer                        footer_type;
	typedef algorithm::CompressionManager       codec_manager_type;
	typedef containers::core::DataBlockSettings settings_type;
	typedef index::Index                        index_type;
	typedef containers::ChecksumContainer       checksum_type;

public:
	TachyonReader();
	TachyonReader(const std::string& filename);
	~TachyonReader();

	/**<
	 * Retrieve current settings records. Settings object is used
	 * prior to each read of a `block` and can be freely modified
	 * before each call. Modifying the settings after a read block
	 * has been invoked has no effect on the loaded `block` data.
	 * @return
	 */
	inline settings_type& getSettings(void){ return(this->settings); }

	/**<
	 * Checks if a FORMAT `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name FORMAT field name to search for (e.g. "GL")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_format_field(const std::string& field_name) const;

	/**<
	 * Checks if a INFO `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name INFO field name to search for (e.g. "AC")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_info_field(const std::string& field_name) const;

	/**<
	 * Checks if a FILTER `field` is set in the header and then checks
	 * if that field exists in the current block. If it does return
	 * the local key. If the field is not described in the header at
	 * all then return -2. If it exists in the header but not in the
	 * loaded block then return -1.
	 * @param field_name FILTER field name to search for (e.g. "PASS")
	 * @return Returns local key if found in this block. Returns -2 if not found in header, or -1 if found in header but not in block
	 */
	const int has_filter_field(const std::string& field_name) const;

	/**<
	 * Calculates which INFO pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name INFO field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_info_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Calculates which FORMAT pattern matches are found for the given field
	 * name in the current loaded block.
	 * @param field_name FORMAT field name
	 * @return           Returns a vector of booleans representing pattern matches
	 */
	const std::vector<bool> get_format_field_pattern_matches(const std::string& field_name) const;

	/**<
	 * Factory function for FORMAT container given an input `field` name
	 * @param field_name FORMAT field name to create a container for
	 * @return           Returns an instance of a `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_format_container(const std::string& field_name) const{
		int format_field = this->has_format_field(field_name);
		if(format_field >= 0) return(new containers::FormatContainer<T>(this->block.format_containers[format_field], this->header.getSampleNumber()));
		else return nullptr;
	}

	/**<
	 * Factory function for balanced FORMAT container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     FORMAT field name to create a container for
	 * @param meta_container Container for meta objects used to balance the FORMAT container
	 * @return               Returns an instance of a balanced `FormatContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::FormatContainer<T>* get_balanced_format_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		int format_field = this->has_format_field(field_name);
		if(format_field >= 0){
			const std::vector<bool> pattern_matches = this->get_format_field_pattern_matches(field_name);
			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			return(new containers::FormatContainer<T>(this->block.format_containers[format_field], meta_container, pattern_matches, this->header.getSampleNumber()));
		}
		else return nullptr;
	}

	/**<
	 * Factory function for INFO container given an input `field` name
	 * @param field_name INFO field name to create a container for
	 * @return           Returns an instance of a `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_info_container(const std::string& field_name) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0) return(new containers::InfoContainer<T>(this->block.info_containers[info_field]));
		else return nullptr;
	}

	/**<
	 * Factory function for balanced INFO container given an input `field` name.
	 * Balances the container such that variants that do not have the given field
	 * will have an empty container placed instead.
	 * @param field_name     INFO field name to create a container for
	 * @param meta_container Container for meta objects used to balance the INFO container
	 * @return               Returns an instance of a balanced `InfoContainer` if successful or a nullpointer otherwise
	 */
	template <class T>
	containers::InfoContainer<T>* get_balanced_info_container(const std::string& field_name, const containers::MetaContainer& meta_container) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0){
			const std::vector<bool> pattern_matches = this->get_info_field_pattern_matches(field_name);

			U32 matches = 0;
			for(U32 i = 0; i < pattern_matches.size(); ++i)
				matches += pattern_matches[i];

			if(matches == 0)
				return nullptr;

			return(new containers::InfoContainer<T>(this->block.info_containers[info_field], meta_container, pattern_matches));
		}
		else return nullptr;
	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void);

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	inline bool open(const std::string& filename){
		this->input_file = filename;
		return(this->open());
	}

	/**<
	 * Overloaded operator for blocks. Useful when
	 * looping over a range of blocks. This happens
	 * frequently in parallel programs.
	 * @param index Block index value in range [0..n_blocks)
	 * @return      Returns TRUE if operation was successful or FALSE otherwise
	 */
	bool operator[](const U32 position);

	/**<
	 *
	 * @param position
	 * @return
	 */
	bool seektoBlock(const U32 position);

	/**<
	 *
	 * @param chromosome_name
	 * @return
	 */
	bool seekToBlockChromosome(const std::string& chromosome_name);

	/**<
	 *
	 * @param chromosome_name
	 * @param from_bp_position
	 * @param to_bp_position
	 * @return
	 */
	bool seekToBlockChromosome(const std::string& chromosome_name, const U32 from_bp_position, const U32 to_bp_position);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool nextBlock(void);

	/**<
	 * Get the current YON block in-order as a copy
	 * @return Returns a YON block. The container has a size of 0 upon fail/empty
	 */
	block_entry_type getBlock(void);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool seek_to_block(const U32& blockID);


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 iterate_all_info(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		// Store as double pointers to avoid memory collisions because
		// infocontainers have different class members
		containers::InfoContainerInterface** its = new containers::InfoContainerInterface*[this->block.index_entry.n_info_streams];
		U32 n_set = 0;
		std::vector<std::string> local_fields;
		for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
			//std::cerr << i << "/" << this->block.index_entry.n_info_streams << ": " << this->block.index_entry.info_offsets[i].key << "@" << this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID << std::endl;
			std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID);
			if(this->header.info_fields[this->block.index_entry.info_offsets[i].key].getType() == core::TACHYON_VARIANT_HEADER_FIELD_TYPE::YON_VCF_HEADER_INTEGER){
				its[n_set++] = new containers::InfoContainer<S32>(this->block.info_containers[i], meta, matches);
				local_fields.push_back(this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID);
			} else if(this->header.info_fields[this->block.index_entry.info_offsets[i].key].getType() == core::TACHYON_VARIANT_HEADER_FIELD_TYPE::YON_VCF_HEADER_STRING ||
                      this->header.info_fields[this->block.index_entry.info_offsets[i].key].getType() == core::TACHYON_VARIANT_HEADER_FIELD_TYPE::YON_VCF_HEADER_CHARACTER){
				its[n_set++] = new containers::InfoContainer<std::string>(this->block.info_containers[i], meta, matches);
				local_fields.push_back(this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID);
			} else if(this->header.info_fields[this->block.index_entry.info_offsets[i].key].getType() == core::TACHYON_VARIANT_HEADER_FIELD_TYPE::YON_VCF_HEADER_FLOAT){
				its[n_set++] = new containers::InfoContainer<double>(this->block.info_containers[i], meta, matches);
				local_fields.push_back(this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID);
			} else {
				//std::cerr << "unknown type: " << (int)this->header.info_fields[this->block.index_entry.info_offsets[i].key].getType() << std::endl;
				its[n_set++] = new containers::InfoContainer<U32>();
				local_fields.push_back(this->header.info_fields[this->block.index_entry.info_offsets[i].key].ID);
				//continue;
			}
		}

		for(U32 p = 0; p < meta.size(); ++p){
			meta.at(p).toVCFString(std::cout, this->header);
			//std::cout << "PASS\t";

			const U32& n_filter_keys = this->block.index_entry.filter_bit_vectors[meta.at(p).filter_pattern_id].n_keys;
			const U32* filter_keys = this->block.index_entry.filter_bit_vectors[meta.at(p).filter_pattern_id].keys;
			if(n_filter_keys){
				std::cout << this->header.filter_fields[filter_keys[0]].ID;
				for(U32 i = 1; i < n_filter_keys; ++i){
					std::cout << ';' << this->header.filter_fields[filter_keys[i]].ID;
				}
				std::cout.put('\t');
			} else
				std::cout << ".\t";

			const U32& n_keys = this->block.index_entry.info_bit_vectors[meta.at(p).info_pattern_id].n_keys;
			const U32* keys = this->block.index_entry.info_bit_vectors[meta.at(p).info_pattern_id].keys;

			// First
			if(this->header.info_fields[this->block.index_entry.info_offsets[keys[0]].key].primitive_type == 2){
				//std::cerr << "skipping" << std::endl;
				std::cout << local_fields[keys[0]];
				continue;
			}
			if(its[keys[0]]->emptyPosition(p)) continue;
			std::cout << local_fields[keys[0]] << "=";
			its[keys[0]]->to_vcf_string(std::cout, p);

			for(U32 i = 1; i < n_keys; ++i){
				//std::cerr << keys[i] << " and " << this->block.index_entry.info_offsets[keys[i]].key << std::endl;
				if(this->header.info_fields[this->block.index_entry.info_offsets[keys[i]].key].primitive_type == 2){
					//std::cerr << "skipping" << std::endl;
					std::cout << ";" << local_fields[keys[i]];
					//exit(1);
					continue;
				}
				if(its[keys[i]]->emptyPosition(p)) continue;
				std::cout << ";" << local_fields[keys[i]] << "=";
				its[keys[i]]->to_vcf_string(std::cout, p);
			}
			std::cout << '\n';
		}

		for(U32 i = 0; i < n_set; ++i)
			delete its[i];

		delete [] its;
		return(0);
	}

	U64 iterate_genotypes(std::ostream& stream = std::cout){
		//std::cerr << "before meta: " << this->block.size() << std::endl;
		containers::MetaContainer meta(this->block);
		//for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i)
		//	std::cerr << this->block.format_containers[i].buffer_data_uncompressed.size() << std::endl;
		//containers::InfoContainer<std::string>* it2 = this->get_balanced_info_container<std::string>("CSQ",meta);
		containers::InfoContainer<U32>* it2 = this->get_balanced_info_container<U32>("DP", meta);
		if(it2 != nullptr){
			std::cerr << it2->size() << std::endl;
		}
		delete it2;

		//std::cerr << meta.size() << std::endl;
		return(meta.size());


	}

	U64 calculateIBS(math::SquareMatrix<double>& square, math::SquareMatrix<double>& square_temporary){
		algorithm::Timer timer;
		timer.Start();

		containers::GenotypeContainer gt(this->block);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].comparePairwise(square_temporary);

		//square /= (U64)2*this->header.getSampleNumber()*gt.size();
		square.addUpperTriagonal(square_temporary, this->block.ppa_manager);
		square_temporary.clear();

		// 2 * (Upper triagonal + diagonal) * number of variants
		const U64 updates = 2*((this->header.getSampleNumber()*this->header.getSampleNumber() - this->header.getSampleNumber())/2 + this->header.getSampleNumber()) * gt.size();
		std::cerr << utility::timestamp("DEBUG") << "Updates: " << utility::ToPrettyString(updates) << '\t' << timer.ElapsedString() << '\t' << utility::ToPrettyString((U64)((double)updates/timer.Elapsed().count())) << "/s" << std::endl;
		return((U64)2*this->header.getSampleNumber()*gt.size());
	}

	U64 getTiTVRatios(std::ostream& stream, std::vector<tachyon::core::TsTvObject>& global){
		containers::GenotypeContainer gt(this->block);
		std::vector<tachyon::core::TsTvObject> objects(this->header.getSampleNumber());
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].getTsTv(objects);

		for(U32 i = 0; i < objects.size(); ++i){
			global[this->block.ppa_manager[i]] += objects[i];
			//stream << i << '\t' << objects[i] << '\n';
		}
		return(gt.size());
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::GenotypeContainer gt(this->block);
		containers::GenotypeSum gt_summary;
		for(U32 i = 0; i < gt.size(); ++i){
			// If there's > 5 alleles continue
			if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
			// Calculate summary statistics
			gt[i].getSummary(gt_summary);

			// Calculate total number of alt-alleles (allele 1, where 0 is ref)
			//std::cerr << gt_summary << '\n';
			gt_summary.clear(); // Recycle summary object
		}
		//std::cerr << std::endl;
		//std::cerr << gt.size() << std::endl;
		return(gt.size());
		//std::cerr << gt[0] << std::endl;;

		//return true;

		core::HeaderMapEntry* entry = nullptr;
		if(this->header.getInfoField("AF", entry)){
			containers::InfoContainer<double> it_i(this->block.info_containers[1]);
			//math::MathSummaryStatistics stats = it_i.getSummaryStatistics();
			//std::cerr << stats.n_total << '\t' << stats.mean << '\t' << stats.standard_deviation << '\t' << stats.min << "-" << stats.max << std::endl;
			for(U32 i = 0; i < it_i.size(); ++i){
				//if(it_i[i].size() < 3) continue;
				//it[i].toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);

				//stream << (int)it_i[i][0];
				for(U32 j = 0; j < it_i[i].size(); ++j)
					stream << it_i[i][j] << ' ';
			}
			stream << '\n';
		}
		return(0);
	}

public:
	std::string        input_file;
	std::ifstream      stream;
	U64                filesize;

	// Actual data
	block_entry_type   block;

	// Supportive objects
	settings_type      settings;
	header_type        header;
	footer_type        footer;
	index_type         index;
	checksum_type      checksums;
	codec_manager_type codec_manager;
};

}

#endif /* CORE_TACHYON_READER_H_ */
