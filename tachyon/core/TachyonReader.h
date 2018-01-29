#ifndef CORE_TACHYONREADER_H_
#define CORE_TACHYONREADER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "../containers/datablock.h"
#include "../algorithm/compression/compression_manager.h"
#include "../algorithm/timer.h"
#include "../index/sorted_index.h"

#include "../containers/abstract_integer_container.h"
#include "../containers/format_container.h"
#include "../containers/genotype_container.h"
#include "../containers/info_container.h"
#include "../containers/meta_cold_container.h"
#include "../containers/meta_hot_container.h"

#include "../iterator/IteratorIntegerReference.h"

#include "../containers/meta_container.h"
#include "../core/genotype_object.h"

#include "../math/fisher.h"
#include "../math/square_matrix.h"

#include "../containers/primitive_group_container.h"
#include "../utility/support_vcf.h"
#include "base/header/yon_tachyonheader.h"

namespace tachyon{

class TachyonReader{
	typedef TachyonReader self_type;
	typedef containers::DataBlock block_entry_type;
	typedef io::BasicBuffer buffer_type;
	typedef core::TachyonHeader header_type;
	typedef algorithm::CompressionManager codec_manager_type;
	typedef containers::core::DataBlockSettings settings_type;
	typedef index::SortedIndex index_type;

public:
	TachyonReader() : filesize(0), n_internal_buffers(0), internal_buffers(nullptr){}
	TachyonReader(const std::string& filename) : input_file(filename), filesize(0), n_internal_buffers(0), internal_buffers(nullptr){}
	~TachyonReader(){
		for(U32 i = 0; i < this->n_internal_buffers; ++i)
			this->internal_buffers[i].deleteAll();

		delete [] this->internal_buffers;
	}

	inline settings_type& getSettings(void){ return(this->settings); }

	/**<
	 * Checks if a target field exists in the header.
	 * @param field_name Field name
	 * @return           Returns the block-offset for that field if it exists. -2 if it does not exist in the header, -1 if it does not exist in the current block
	 */
	inline const int has_format_field(const std::string& field_name) const{
		core::HeaderMapEntry* match = nullptr;
		if(this->header.getEntry(field_name, match)){
			U32 target = -1;
			for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i){
				//std::cerr << i << '/' << this->block.index_entry.n_format_streams << '\t' << this->block.index_entry.format_offsets[i].key << '\t' << this->header.entries[this->block.index_entry.format_offsets[i].key].ID << std::endl;
				if(this->block.index_entry.format_offsets[i].key == match->IDX){
					target = i;
					break;
				}
			}
			//std::cerr << "target stream is: " << target << std::endl;
			return(target);
		}
		return(-2);
	}

	inline const int has_info_field(const std::string& field_name) const{
		core::HeaderMapEntry* match = nullptr;
		if(this->header.getEntry(field_name, match)){
			U32 target = -1;
			for(U32 i = 0; i < this->block.index_entry.n_info_streams; ++i){
				//std::cerr << i << '/' << this->block.index_entry.n_info_streams << '\t' << this->block.index_entry.info_offsets[i].key << '\t' << this->header.entries[this->block.index_entry.info_offsets[i].key].ID << std::endl;
				if(this->block.index_entry.info_offsets[i].key == match->IDX){
					target = i;
					break;
				}
			}
			//std::cerr << "target stream is: " << target << std::endl;
			return(target);
		}
		return(-2);
	}

	inline const std::vector<bool> get_info_field_pattern_matches(const std::string& field_name) const{
		int info_field = this->has_info_field(field_name);
		std::vector<bool> ret;
		if(info_field >= 0){
			// Collect all matches
			// Place in array
			// 0 = false, 1 = true
			ret.resize(this->block.index_entry.n_info_patterns, false);
			for(U32 i = 0; i < this->block.index_entry.n_info_patterns; ++i){
				std::cerr << i << '\t' << this->block.index_entry.info_bit_vectors[i][info_field] << std::endl;
				ret[i] = this->block.index_entry.info_bit_vectors[i][info_field];
			}
		}
		return(ret);
	}

	template <class T>
	containers::FormatContainer<T>* get_format_container(const std::string& field_name) const{
		int format_field = this->has_format_field(field_name);
		if(format_field >= 0) return(new containers::FormatContainer<T>(this->block.format_containers[format_field], this->header.n_samples));
		else return nullptr;
	}

	template <class T>
	containers::InfoContainer<T>* get_info_container(const std::string& field_name) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0) return(new containers::InfoContainer<T>(this->block.info_containers[info_field]));
		else return nullptr;
	}

	// Todo: balance container with 0's when data is missing
	template <class T>
	containers::InfoContainer<T>* get_balanced_info_container(const std::string& field_name) const{
		int info_field = this->has_info_field(field_name);
		if(info_field >= 0) return(new containers::InfoContainer<T>(this->block.info_containers[info_field]));
		else return nullptr;
	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(void){
		if(this->input_file.size() == 0){
			std::cerr << "no filename" << std::endl;
			return false;
		}
		this->stream.open(this->input_file, std::ios::binary | std::ios::in | std::ios::ate);
		this->filesize = (U64)this->stream.tellg();
		if(!this->stream.good()){
			std::cerr << "failed to read file" << std::endl;
			return false;
		}
		this->stream.seekg(0);
		if(!this->stream.good()){
			std::cerr << "failed to rewind" << std::endl;
			return false;
		}

		// Load header
		this->stream << this->header;
		if(!this->stream.good()){
			std::cerr << "failed to get header" << std::endl;
			return false;
		}
		// Keep track of start position
		const U64 return_pos = this->stream.tellg();

		// Seek to EOF and make check
		this->stream.seekg(this->filesize - 32);
		BYTE eof_data[32];
		helpers::HexToBytes(constants::TACHYON_FILE_EOF, &eof_data[0]);
		BYTE eof_match[32];
		this->stream.read((char*)&eof_match[0], 32);
		for(U32 i = 0; i < 32; ++i){
			if(eof_data[i] != eof_match[i]){
				std::cerr << "File is truncated!" << std::endl;
				return false;
			}
		}
		// Seek back to find start of index
		// Seek to that start of index
		// Load index
		// Seek back to start of the file
		this->stream.seekg(this->filesize - 32 - sizeof(U64));
		this->stream.read((char*)reinterpret_cast<char*>(&this->l_data), sizeof(U64));
		this->stream.seekg(this->l_data);
		this->stream >> this->index;
		this->stream.seekg(return_pos);

		return(this->stream.good());
	}

	/**<
	 * Opens a YON file. Performs all prerequisite
	 * checks and loads all auxiliary data structures
	 * @param filename Target input filename
	 * @return Returns TRUE upon success or FALSE otherwise
	 */
	bool open(const std::string& filename){
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
	bool operator[](const U32 index);

	/**<
	 * Get the next YON block in-order
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool get_next_block(){
		// If the stream is faulty then return
		if(!this->stream.good()){
			std::cerr << "faulty stream" << std::endl;
			return false;
		}

		// If the current position is the EOF then
		// exit the function
		if((U64)this->stream.tellg() == this->l_data)
			return false;

		// Reset and re-use
		this->block.clear();

		// Attempts to read a YON block with the provided
		// settings
		if(!this->block.read(stream, this->settings))
			return false;

		// Internally decompress available data
		if(!this->codec_manager.decompress(this->block))
			return false;

		// All passed
		return true;
	}

	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param b
	 * @return
	 */
	bool seek_to_block(const U32& b);

	U64 iterate_genotypes(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		std::cerr << block.size() << std::endl;

		std::vector<bool> ac_matches = this->get_info_field_pattern_matches("AC");
		return(0);

		containers::InfoContainer<U32>* af = this->get_info_container<U32>("AC");
		if(af != nullptr){
			std::cerr << "in info container: " << af->size() << std::endl;
			for(U32 variant = 0; variant < af->size(); ++variant){
				util::to_vcf_string(stream, (*af)[variant]) << ' ';
			}
			std::cerr << std::endl;
		} else std::cerr << "AC not found" << std::endl;
		delete af;
		return(0);

		containers::FormatContainer<float>* it = this->get_format_container<float>("GP");
		if(it != nullptr){
			for(U32 variant = 0; variant < it->size(); ++variant){ // variants
				for(U32 sample = 0; sample < (*it)[variant].size(); ++sample){ // individuals
					util::to_vcf_string(std::cout, (*it)[variant][sample]);
					std::cerr<<"\t";
				}
				std::cerr << '\n';
			}
			std::cerr << std::endl;
		}
		delete it;
		return this->block.size();
	}

	void calculateIBS(math::SquareMatrix<double>& square){
		algorithm::Timer timer;
		timer.Start();

		containers::GenotypeContainer gt(this->block);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].compareSamplesPairwise(square);

		square /= (U64)2*this->header.n_samples*gt.size();
		const U64 updates = 2*(this->header.n_samples*this->header.n_samples - this->header.n_samples)/2 * gt.size();
		std::cerr << "Updates: " << helpers::ToPrettyString(updates) << '\t' << timer.ElapsedString() << '\t' << helpers::ToPrettyString((U64)((double)updates/timer.Elapsed().count())) << "/s" << std::endl;
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::GenotypeContainer gt(this->block);
		math::Fisher fisher(1000);
		containers::GenotypeSum gt_summary;
		for(U32 i = 0; i < gt.size(); ++i){
			//std::vector<core::GTObject> objects = gt[i].getObjects();
			//const U32 n_entries = gt[i].getSum();
			if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
			gt[i].getSummary(gt_summary);

			const U64 total = gt_summary.getAlleleA(1) + gt_summary.getAlleleB(1);
			const double p = fisher.fisherTest(gt_summary.getAlleleA(1), total, gt_summary.getAlleleB(1), total);
			//if(p < 1e-3){
			//	gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
			//	stream << '\t' << gt_summary << '\t' << p << '\t' << ((gt_summary.getAlleleA(1) == 0 || gt_summary.getAlleleB(1) == 0) ? 1 : 0) << '\n';
			//}

			//if(!gt[i].getMeta().isBiallelic()){
			//	gt[i].getMeta().toVCFString(stream, this->header, this->block.index_entry.contigID, this->block.index_entry.minPosition);
			//	std::cerr << '\t' << gt_summary << std::endl;
			//}
			//assert(gt_summary.alleleCount() == 2*this->header.n_samples);
			//assert(gt_summary.genotypeCount() == this->header.n_samples);
			//assert(n_entries == this->header.n_samples);
			gt_summary.clear();
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

#endif /* CORE_TACHYONREADER_H_ */
