#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/encryption/EncryptionDecorator.h"
#include "algorithm/timer.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/primitive_group_container.h"
#include "containers/meta_container.h"
#include "algorithm/digital_digest.h"
#include "containers/variantblock.h"
#include "core/genotype_object.h"
#include "core/footer/footer.h"
#include "core/header/variant_header.h"
#include "math/fisher.h"
#include "math/square_matrix.h"
#include "math/basic_vector_math.h"
#include "utility/support_vcf.h"
#include "index/index.h"

namespace tachyon{

struct VariantReaderObjects{
public:
	typedef VariantReaderObjects                 self_type;
	typedef containers::MetaContainer            meta_container_type;
	typedef containers::GenotypeContainer        gt_container_type;
	typedef containers::InfoContainerInterface   info_interface_type;
	typedef containers::FormatContainerInterface format_interface_type;

public:
	VariantReaderObjects() :
		loaded_genotypes(false),
		loaded_meta(false),
		n_loaded_info(0),
		n_loaded_format(0),
		meta(nullptr),
		genotypes(nullptr),
		info_fields(nullptr),
		format_fields(nullptr)
	{}

	~VariantReaderObjects(){
		delete this->meta;
		delete this->genotypes;

		for(U32 i = 0; i < this->n_loaded_info; ++i)
			delete this->info_fields[i];
		delete [] this->info_fields;

		for(U32 i = 0; i < this->n_loaded_format; ++i)
			delete this->format_fields[i];
		delete [] this->format_fields;
	}

public:
	bool loaded_genotypes;
	bool loaded_meta;
	size_t n_loaded_info;
	size_t n_loaded_format;
	std::vector<U32> info_keep;
	std::vector<U32> format_keep;
	std::vector< std::vector<U32> > local_match_keychain_info;
	std::vector< std::vector<U32> > local_match_keychain_format;

	std::vector<std::string> info_field_names;
	std::vector<std::string> format_field_names;

	meta_container_type*     meta;
	gt_container_type*       genotypes;
	info_interface_type**    info_fields;
	format_interface_type**  format_fields;
};

class VariantReader{
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef core::VariantHeader                    header_type;
	typedef core::Footer                           footer_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef core::DataBlockSettings                settings_type;
	typedef index::Index                           index_type;
	typedef algorithm::VariantDigitalDigestManager checksum_type;
	typedef encryption::Keychain                   keychain_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;

	// Function pointers
	typedef void (self_type::*print_format_function)(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	typedef void (self_type::*print_info_function)(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	typedef void (self_type::*print_filter_function)(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	typedef buffer_type& (*print_meta_function)(buffer_type& buffer, const char& delimiter, const meta_entry_type& meta_entry, const header_type& header, const core::SettingsCustomOutput& controller);

public:
	VariantReader();
	VariantReader(const std::string& filename);
	VariantReader(const self_type& other);
	~VariantReader();

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

	/**<
	 *
	 * @return
	 */
	bool parseSettings(void){
		settings.load_info_ID_loaded.clear();
		settings.load_format_ID_loaded.clear();

		// Map INFO
		if(settings.load_info == false){ // prevent double load
			U32 info_matches = 0;
			for(U32 i = 0; i < settings.info_list.size(); ++i){
				const S32 global_key = this->has_info_field(settings.info_list[i]);
				if(global_key >= 0){
					S32 local_key = -1;
					for(U32 i = 0; i < this->block.footer.n_info_streams; ++i){
						if(this->block.footer.info_offsets[i].data_header.global_key == global_key){
							local_key = i;
							this->settings.info_ID_list.push_back(i);
							settings.load_info_ID_loaded.push_back(
														core::SettingsMap(
																info_matches++, // iterator value
																i,              // local index id
																&this->block.footer.info_offsets[i]) // offset
																);
							break;
						}
					}

					if(local_key == -1){
						//std::cerr << "could not find local" << std::endl;
					}
				}
			}
		}

		// Map FORMAT
		if(settings.load_format == false){ // prevent double load
			U32 format_matches = 0;
			for(U32 i = 0; i < settings.format_list.size(); ++i){
				const S32 global_key = this->has_format_field(settings.format_list[i]);
				if(global_key >= 0){
					S32 local_key = -1;
					for(U32 i = 0; i < this->block.footer.n_format_streams; ++i){
						if(this->block.footer.format_offsets[i].data_header.global_key == global_key){
							local_key = i;
							this->settings.format_ID_list.push_back(i);
							settings.load_format_ID_loaded.push_back(
														core::SettingsMap(
																format_matches++, // iterator value
																i,              // local index id
																&this->block.footer.format_offsets[i]) // offset
																);
							break;
						}
					}

					if(local_key == -1){
						//std::cerr << "could not find local" << std::endl;
					}
				}
			}
		}

		return(true);
	}

	/**<
	 * Primary construction function for generating the appropriate instances of
	 * iterators / containers
	 * @param objects Target objects
	 * @return        Returns reference to input target objects
	 */
	objects_type& loadObjects(objects_type& objects) const;

	/**<
	 * Outputs
	 * @return
	 */
	const U64 outputVCF(void){
		U64 n_variants = 0;

		// Output VCF header
		if(this->settings.show_vcf_header){
			this->header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
			this->header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
					  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

			this->header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE;

			this->header.writeVCFHeaderString(std::cout, this->settings.load_format || this->settings.format_list.size());
		}

		// While there are YON blocks
		while(this->nextBlock()) n_variants += this->outputBlockVCF();
		return(n_variants);
	}

	/**<
	 *
	 * @return
	 */
	const U64 outputCustom(void){
		U64 n_variants = 0;

		// While there are YON blocks
		while(this->nextBlock()) n_variants += this->outputBlockCustom();
		return(n_variants);
	}

	/**<
	 *
	 * @return
	 */
	const U32 outputBlockVCF(void) const{
		objects_type objects;
		this->loadObjects(objects);

		// Reserve memory for output buffer
		// This is much faster than writing directly to ostream because of syncing
		io::BasicBuffer output_buffer(256000 + this->header.getSampleNumber()*2);
		std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());

		for(U32 p = 0; p < objects.meta->size(); ++p){
			if(this->settings.custom_output_format)
				utility::to_vcf_string(output_buffer, '\t', (*objects.meta)[p], this->header, this->settings.custom_output_controller);
			else
				utility::to_vcf_string(output_buffer, '\t', (*objects.meta)[p], this->header);

			// Filter options
			if(settings.load_set_membership) this->printFILTER(output_buffer, p, objects);
			else output_buffer += '.';

			if(settings.load_info || settings.load_format || this->block.n_info_loaded) output_buffer += '\t';
			else {
				output_buffer += '\n';
				continue;
			}

			// Print normal or with a custom delimiter
			this->printINFOVCF(output_buffer, p, objects);

			if(settings.load_format || this->block.n_format_loaded){
				output_buffer += this->settings.custom_delimiter_char;

				if(this->settings.format_ID_list.size())
					this->printFORMATCustom(output_buffer, '\t', p, objects, genotypes_unpermuted);
				else
					this->printFORMATVCF(output_buffer, p, objects, genotypes_unpermuted);

			}

			output_buffer += '\n';

			// Add FORMAT data

			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

		return(objects.meta->size());
	}

	/**<
	 *
	 * @return
	 */
	const U32 outputBlockCustom(void) const{
		objects_type objects;
		this->loadObjects(objects);

		// Reserve memory for output buffer
		// This is much faster than writing directly to ostream because of syncing
		io::BasicBuffer output_buffer(256000 + this->header.getSampleNumber()*2);
		std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());

		// Todo: move to function
		U32 info_match_limit = 1; // any match
		//info_match_limit = this->settings.info_list.size(); // all match

		// Function pointer to use
		print_format_function print_format = &self_type::printFORMATDummy;
		print_info_function   print_info   = &self_type::printINFODummy;
		print_meta_function   print_meta   = &utility::to_vcf_string;
		print_filter_function print_filter = &self_type::printFILTERDummy;

		if(settings.output_json) print_meta = &utility::to_json_string;

		if(settings.load_format || this->block.n_format_loaded){
			if(settings.output_json){
				print_format = &self_type::printFORMATCustomVectorJSON;
			} else {
				if(settings.output_format_vector) print_format = &self_type::printFORMATCustomVector;
				else print_format = &self_type::printFORMATCustom;
			}
		}

		if(settings.load_info || this->block.n_info_loaded){
			if(settings.output_json) print_info = &self_type::printINFOCustomJSON;
			else print_info = &self_type::printINFOCustom;
		}

		if(settings.custom_output_controller.show_filter){
			if(settings.output_json) print_filter = &self_type::printFILTERJSON;
			else print_filter = &self_type::printFILTERCustom;
		}

		U32 n_records_returned = 0;

		if(settings.output_json) output_buffer += "\"block\":{";
		for(U32 position = 0; position < objects.meta->size(); ++position){
			//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
			//	continue;
			if(settings.output_json){
				if(position != 0) output_buffer += ",\n";

				output_buffer += "\"obj-";
				output_buffer.AddReadble(position);
				output_buffer += "\":{";
			}
			++n_records_returned;

			(*print_meta)(output_buffer, this->settings.custom_delimiter_char, (*objects.meta)[position], this->header, this->settings.custom_output_controller);
			(this->*print_filter)(output_buffer, position, objects);
			(this->*print_info)(output_buffer, this->settings.custom_delimiter_char, position, objects);
			(this->*print_format)(output_buffer, this->settings.custom_delimiter_char, position, objects, genotypes_unpermuted);

			if(settings.output_json) output_buffer += "}";
			else output_buffer += '\n';
			//output_buffer += "}";

			// Flush if buffer is large
			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}
		if(settings.output_json) output_buffer += "}";

		// Flush buffer
		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

		exit(1);
		return(n_records_returned);
	}

	// Dummy functions as interfaces for function pointers
	inline void printFILTERDummy(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const{}
	inline void printFORMATDummy(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const{}
	inline void printINFODummy(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const{}

	// FILTER functions
	void printFILTER(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printFILTERCustom(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printFILTERJSON(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;

	// FORMAT functions
	void printFORMATVCF(buffer_type& buffer, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATVCF(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustom(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustomVector(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	void printFORMATCustomVectorJSON(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;


	// INFO functions
	void printINFOVCF(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	void printINFOVCF(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	void printINFOCustom(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	void printINFOCustomJSON(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 timings_meta(){
		containers::MetaContainer meta(this->block);
		buffer_type temp(meta.size() * 1000);

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(temp, this->settings.custom_delimiter_char, meta[p], this->header);
			//utility::to_vcf_string(std::cout, meta[p], this->header);
			temp += '\n';
			//if(temp.size() > 65536){
			//	std::cout.write(temp.data(), temp.size());
			//	temp.reset();
			//}
		}
		std::cout.write(temp.data(), temp.size());
		return(meta.size());
	}


	U64 iterate_genotypes(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);

		for(U32 i = 0; i < gt.size(); ++i){
			// All of these functions are in relative terms very expensive!
			// Avoid using them unless you absolutely have to!
			// Vector of literal genotype representations (lower level)
			//std::vector<core::GTObject> objects     = gt[i].getLiteralObjects();
			// Vector of genotype objects (high level permuted)
			//std::vector<core::GTObject> objects_all = gt[i].getObjects(this->header.getSampleNumber());
			// Vector of genotype objects (high level unpermuted - original)
			std::vector<core::GTObject> objects_true = gt[i].getObjects(this->header.getSampleNumber(), this->block.ppa_manager);

			std::cout << (int)objects_true[i].alleles[0].first << (objects_true[i].alleles[1].second ? '/' : '|') << (int)objects_true[i].alleles[1].first;
			for(U32 i = 1; i < objects_true.size(); ++i){
				std::cout << '\t' << (int)objects_true[i].alleles[0].first << (objects_true[i].alleles[1].second ? '/' : '|') << (int)objects_true[i].alleles[1].first;
			}
			std::cout << std::endl;
		}
		return(gt.size());
	}

	U64 calculateIBS(math::SquareMatrix<double>& square, math::SquareMatrix<double>& square_temporary){
		algorithm::Timer timer;
		timer.Start();

		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);
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

	U64 getTiTVRatios(std::ostream& stream, std::vector<core::TsTvObject>& global){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);

		std::vector<core::TsTvObject> objects(this->header.getSampleNumber());
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].getTsTv(objects);

		for(U32 i = 0; i < objects.size(); ++i)
			global[this->block.ppa_manager[i]] += objects[i];

		return(gt.size());
	}

	U64 getGenotypeSummary(std::ostream& stream){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);
		containers::GenotypeSummary gtsum(10);

		for(U32 i = 0; i < gt.size(); ++i){
			gt[i].getSummary(gtsum);
			std::vector<double> hwe_p = gtsum.calculateHardyWeinberg(meta[i]);
			std::vector<double> af = gtsum.calculateAlleleFrequency(meta[i]);
			//if(hwe_p[0] < 1e-3){
				utility::to_vcf_string(stream, this->settings.custom_delimiter_char, meta[i], this->header);
				stream << "AF=" << af[0];
				for(U32 p = 1; p < af.size(); ++p){
					stream << "," << af[p];
				}
				stream << ";HWE_P=" << hwe_p[0];

				for(U32 p = 1; p < hwe_p.size(); ++p){
					stream << "," << hwe_p[p];
				}
				stream.put('\n');
			//}
			gtsum.clear();
		}
		return(gt.size());
	}

	U64 countVariants(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		return(meta.size());
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer gt(this->block, meta);
		containers::GenotypeSummary gt_summary;
		for(U32 i = 0; i < gt.size(); ++i){
			// If there's > 5 alleles continue
			if(gt[i].getMeta().getNumberAlleles() >= 5) continue;
			// Calculate summary statistics
			//gt[i].getSummary(gt_summary);

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
	keychain_type      keychain;
};

}

#endif /* CORE_TACHYON_READER_H_ */
