#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/timer.h"
#include "containers/datablock.h"
#include "containers/abstract_integer_container.h"
#include "containers/format_container.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/meta_cold_container.h"
#include "containers/meta_hot_container.h"
#include "containers/primitive_group_container.h"
#include "containers/meta_container.h"
#include "core/genotype_object.h"
#include "core/header/tachyon_header.h"
#include "math/fisher.h"
#include "math/square_matrix.h"
#include "math/basic_vector_math.h"
#include "utility/support_vcf.h"
#include "iterator/IteratorIntegerReference.h"
#include "index/sorted_index.h"

namespace tachyon{

class TachyonReader{
	typedef TachyonReader                       self_type;
	typedef containers::DataBlock               block_entry_type;
	typedef io::BasicBuffer                     buffer_type;
	typedef core::TachyonHeader                 header_type;
	typedef algorithm::CompressionManager       codec_manager_type;
	typedef containers::core::DataBlockSettings settings_type;
	typedef index::SortedIndex                  index_type;

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
		if(format_field >= 0) return(new containers::FormatContainer<T>(this->block.format_containers[format_field], this->header.n_samples));
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

			return(new containers::FormatContainer<T>(this->block.format_containers[format_field], meta_container, pattern_matches, this->header.n_samples));
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
	bool get_next_block();

	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param b
	 * @return
	 */
	bool seek_to_block(const U32& position);

	//<----------------- EXAMPLE FUNCTIONS -------------------------->
	U64 iterate_genotypes(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		//std::cerr << block.size() << std::endl;

		// Variant-balanced
		containers::InfoContainer<char>* it2 = this->get_balanced_info_container<char>("CS", meta);
		// Not variant-balanced
		//containers::InfoContainer<U32>* it3 = this->get_info_container<U32>("AR2");

		//containers::FormatContainer<float>* it4 = this->get_balanced_format_container<float>("GP", meta);
		if(it2 != nullptr){
			//std::cerr << "balanced format = " << it4->size() << std::endl;
			//math::SummaryStatistics ss = math::summary_statistics(*it4);
			//std::cerr << ss.mean << "+-" << ss.getSigmaSquared() << " (" << ss.min << "-" << ss.max << ") n = " << ss.n_total << std::endl;

			for(size_t i = 0; i < it2->size(); ++i){
				if(it2->at(i).size() == 0) continue;
				else {

					meta.at(i).toVCFString(std::cout, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
					std::cout << "\t";
					utility::to_vcf_string(std::cout, (*it2)[i]);
					std::cout << "\n";
				}


			}

			/*
			for(U32 i = 0; i < it4->size(); ++i){
				for(U32 j = 0; j < it4->at(i).size(); ++j){
					//util::to_vcf_string(stream, it4->at(i).at(j)) << ' ';
					//std::cerr << math::max(it4->at(i).at(j)) << ' ';
					const math::SummaryStatistics ss = math::summary_statistics(it4->at(i).at(j));
					stream << ss.min << "-" << ss.max << "(" << ss.mean << "," << ss.getSigmaSquared() << ") ";
				}
				stream << '\n';
			}
			stream << '\n';
			*/
		}
		delete it2;

		//if(it2!=nullptr)   std::cerr << "it  = " << it2->size() << std::endl;
		//if(it3 != nullptr) std::cerr << "it2 = " << it3->size() << std::endl;
		//delete it2;
		return(0);

		containers::InfoContainer<U32>* af = this->get_info_container<U32>("AC");
		if(af != nullptr){
			std::cerr << "in info container: " << af->size() << std::endl;
			for(U32 variant = 0; variant < af->size(); ++variant){
				utility::to_vcf_string(stream, (*af)[variant]) << ' ';
			}
			std::cerr << std::endl;
		} else std::cerr << "AC not found" << std::endl;
		delete af;
		return(0);

		containers::FormatContainer<float>* it = this->get_format_container<float>("GP");
		if(it != nullptr){
			for(U32 variant = 0; variant < it->size(); ++variant){ // variants
				for(U32 sample = 0; sample < (*it)[variant].size(); ++sample){ // individuals
					utility::to_vcf_string(std::cout, (*it)[variant][sample]);
					std::cerr<<"\t";
				}
				std::cerr << '\n';
			}
			std::cerr << std::endl;
		}
		delete it;
		return this->block.size();
	}

	U64 calculateIBS(math::SquareMatrix<double>& square, math::SquareMatrix<double>& square_temporary){
		algorithm::Timer timer;
		timer.Start();

		containers::GenotypeContainer gt(this->block);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].comparePairwise(square_temporary);

		//square /= (U64)2*this->header.n_samples*gt.size();
		square.addUpperTriagonal(square_temporary, this->block.ppa_manager);
		square_temporary.clear();

		// 2 * (Upper triagonal + diagonal) * number of variants
		const U64 updates = 2*((this->header.n_samples*this->header.n_samples - this->header.n_samples)/2 + this->header.n_samples) * gt.size();
		std::cerr << utility::timestamp("DEBUG") << "Updates: " << utility::ToPrettyString(updates) << '\t' << timer.ElapsedString() << '\t' << utility::ToPrettyString((U64)((double)updates/timer.Elapsed().count())) << "/s" << std::endl;
		return((U64)2*this->header.n_samples*gt.size());
	}

	U64 getTiTVRatios(std::ostream& stream, std::vector<tachyon::core::TiTvObject>& global){
		containers::GenotypeContainer gt(this->block);
		std::vector<tachyon::core::TiTvObject> objects(this->header.n_samples);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].updateTransitionTransversions(objects);

		for(U32 i = 0; i < objects.size(); ++i){
			global[this->block.ppa_manager[i]] += objects[i];
			//stream << i << '\t' << objects[i] << '\n';
		}
		return(gt.size());
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::GenotypeContainer gt(this->block);
		//math::Fisher fisher(1000);
		//containers::GenotypeSum gt_summary1, gt_summary2;
		for(U32 i = 0; i < gt.size(); ++i){
			//std::vector<core::GTObject> objects = gt[i].getLiteralObjects();
			//std::vector<core::GTObject> objects_all  = gt[i].getObjects(this->header.n_samples);
			std::vector<core::GTObject> objects_true = gt[i].getObjects(this->header.n_samples, this->block.ppa_manager);
			//std::cerr << objects.size() << '\t' << objects_all.size() << std::endl;
			gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
			//std::cerr << "\nPermuted: ";
			//for(U32 j = 0; j < objects_all.size(); ++j)
			//	std::cerr << (int)objects_all[j].alleles[0].first << (objects_all[j].alleles[0].second ? "|" : "/") << (int)objects_all[j].alleles[1].first << ' ';
			//std::cerr << std::endl;

			/*
			algorithm::PermutationManager& ppa = this->block.ppa_manager;
			core::GTObject** pointers = new core::GTObject*[objects_all.size()];
			//for(U32 j = 0; j < objects_all.size(); ++j) pointers[j] = nullptr;
			for(U32 j = 0; j < objects_all.size(); ++j){
				//assert(pointers[ppa[j]] == nullptr);
				pointers[ppa[j]] = &objects_all[j];
			}

			//for(U32 j = 0; j < objects_all.size(); ++j) assert(pointers[j] != nullptr);

			//std::cerr << "Restored: ";
			for(U32 j = 0; j < objects_all.size(); ++j){
				stream << (int)pointers[j]->alleles[0].first << (pointers[j]->alleles[1].second ? "|" : "/") << (int)pointers[j]->alleles[1].first << ' ';
			}
			stream << '\n';

			delete [] pointers;
			//delete pointers;
	*/

			utility::to_vcf_string(stream, objects_true) << '\n';
			//utility::to_vcf_string(stream, objects_true[0]) << '\t';
			/*
			for(U32 j = 1; j < objects_true.size(); ++j){
				stream << '\t';
				utility::to_vcf_string(stream, objects_true[j]);
				//stream << (int)objects_true[j].alleles[0].first << (objects_true[j].alleles[1].second ? "|" : "/") << (int)objects_true[j].alleles[1].first << ' ';
			}
			stream << '\n';
			*/
			/*

			const U32 n_entries = gt[i].getSum();
			if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
			gt[i].getSummary(gt_summary);

			const U64 total = gt_summary.getAlleleA(1) + gt_summary.getAlleleB(1);
			const double p = fisher.fisherTest(gt_summary.getAlleleA(1), total, gt_summary.getAlleleB(1), total);
			if(p < 1e-3){
				gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
				stream << '\t' << gt_summary << '\t' << p << '\t' << ((gt_summary.getAlleleA(1) == 0 || gt_summary.getAlleleB(1) == 0) ? 1 : 0) << '\n';
			}

			//if(!gt[i].getMeta().isBiallelic()){
			//	gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
			//	std::cerr << '\t' << gt_summary << std::endl;
			//}
			assert(gt_summary.alleleCount() == 2*this->header.n_samples);
			assert(gt_summary.genotypeCount() == this->header.n_samples);
			assert(n_entries == this->header.n_samples);
			gt_summary.clear();
			*/
		}
		//std::cerr << std::endl;
		//std::cerr << gt.size() << std::endl;
		return(gt.size());
		//std::cerr << gt[0] << std::endl;;

		//return true;

		core::HeaderMapEntry* entry = nullptr;
		if(this->header.getEntry("AF", entry)){
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
	std::string   input_file;
	std::ifstream stream;
	U64           filesize;
	U64           l_data;   // Size in bytes from end of header to start of footer

	settings_type      settings;
	header_type        header;
	index_type         index;
	block_entry_type   block;
	codec_manager_type codec_manager;

private:
	U32          n_internal_buffers;
	buffer_type* internal_buffers;
};

}

#endif /* CORE_TACHYON_READER_H_ */
