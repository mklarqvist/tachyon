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
	objects_type& loadObjects(objects_type& objects) const{
		objects.meta = new meta_container_type(this->block);
		if(this->block.header.controller.hasGT && settings.load_genotypes_all){
			objects.genotypes = new gt_container_type(this->block, *objects.meta);
		}

		objects.n_loaded_format = this->block.n_format_loaded;
		objects.format_fields = new format_interface_type*[objects.n_loaded_format];

		if(this->block.n_format_loaded){
			for(U32 i = 0; i < this->block.n_format_loaded; ++i){
				const U32 global_key = settings.load_format_ID_loaded[i].offset->data_header.global_key;
				std::vector<bool> matches = this->get_format_field_pattern_matches(this->header.format_fields[global_key].ID);

				if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					objects.format_fields[i] = new containers::FormatContainer<S32>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
					objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					objects.format_fields[i] = new containers::FormatContainer<std::string>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
					objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					objects.format_fields[i] = new containers::FormatContainer<float>(this->block.format_containers[i], *objects.meta, matches, this->header.getSampleNumber());
					objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
				} else {
					objects.format_fields[i] = new containers::FormatContainer<U32>;
					objects.format_field_names.push_back(this->header.format_fields[global_key].ID);
				}
			}
		}

		// Store as double pointers to avoid memory collisions because
		// info containers have different class members
		objects.n_loaded_info = this->block.n_info_loaded;
		objects.info_fields = new info_interface_type*[objects.n_loaded_info];

		if(this->block.n_info_loaded){
			for(U32 i = 0; i < this->block.n_info_loaded; ++i){
				const U32 global_key = settings.load_info_ID_loaded[i].offset->data_header.global_key;
				std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[global_key].ID);

				if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					objects.info_fields[i] = new containers::InfoContainer<S32>(this->block.info_containers[i], *objects.meta, matches);
					objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					objects.info_fields[i] = new containers::InfoContainer<std::string>(this->block.info_containers[i], *objects.meta, matches);
					objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					objects.info_fields[i] = new containers::InfoContainer<float>(this->block.info_containers[i], *objects.meta, matches);
					objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
				} else {
					objects.info_fields[i] = new containers::InfoContainer<U32>();
					objects.info_field_names.push_back(this->header.info_fields[global_key].ID);
				}
			}
		}

		// If we want to drop records that do not have all/any of the fields we desire
		// then we create a vector of size N_PATTERNS and set those that MATCH to TRUE
		// this allows for filtering in O(1)-time
		//
		// This vector stores the number of INFO fields having set membership with this
		// particular hash pattern
		objects.info_keep = std::vector<U32>(this->block.footer.n_info_patterns, 0);
		objects.format_keep = std::vector<U32>(this->block.footer.n_format_patterns, 0);

		// This vector of vectors keeps the local INFO identifiers for the matched global
		// identifiers in a given hash pattern
		//
		// For example: x[5] contains the local IDs for loaded INFO streams for pattern ID 5
		objects.local_match_keychain_info = std::vector< std::vector<U32> >(this->block.footer.n_info_patterns);
		objects.local_match_keychain_format = std::vector< std::vector<U32> >(this->block.footer.n_format_patterns);

		// If loading all INFO values then return them in the ORIGINAL order
		if(this->settings.load_info){
			for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // Number of info patterns
				for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
					for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
						if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
							objects.local_match_keychain_info[i].push_back(k);
							++objects.info_keep[i];
						}
					}
				}
			}
		}
		// If loading custom INFO fields then return them in the REQUESTED order
		else {
			for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // i = Number of info patterns
				for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
					for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
						if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
							objects.local_match_keychain_info[i].push_back(k);
							++objects.info_keep[i];
						}
					}
				}
			}
		}

		// For FORMAT
		// If loading all FORMAT values then return them in the ORIGINAL order
		if(this->settings.load_format){
			for(U32 i = 0; i < this->block.footer.n_format_patterns; ++i){ // Number of info patterns
				for(U32 j = 0; j < this->block.footer.format_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
					for(U32 k = 0; k < settings.load_format_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
						if(this->block.footer.format_offsets[this->block.footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_format_ID_loaded[k].offset->data_header.global_key){
							objects.local_match_keychain_format[i].push_back(k);
							++objects.format_keep[i];
						}
					}
				}
			}
		}
		// If loading custom FORMAT fields then return them in the REQUESTED order
		else {
			for(U32 i = 0; i < this->block.footer.n_format_patterns; ++i){ // i = Number of info patterns
				for(U32 k = 0; k < settings.load_format_ID_loaded.size(); ++k){ // k = Number of loaded INFO identifiers
					for(U32 j = 0; j < this->block.footer.format_bit_vectors[i].n_keys; ++j){ // j = Number of keys in pattern [i]
						if(this->block.footer.format_offsets[this->block.footer.format_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_format_ID_loaded[k].offset->data_header.global_key){
							objects.local_match_keychain_format[i].push_back(k);
							++objects.format_keep[i];
						}
					}
				}
			}
		}

		return(objects);
	}

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

			this->header.writeVCFHeaderString(std::cout, this->settings.load_format);
		}

		// While there are YON blocks
		while(this->nextBlock()){
			n_variants += this->outputBlockVCF();
		}
		return(n_variants);
	}

	/**<
	 *
	 * @return
	 */
	const U64 outputCustom(void){
		U64 n_variants = 0;

		// While there are YON blocks
		while(this->nextBlock()){
			n_variants += this->outputBlockCustom();
		}
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
				utility::to_vcf_string(output_buffer, this->settings.custom_delimiter_char, (*objects.meta)[p], this->header, this->settings.custom_output_controller);
			else
				utility::to_vcf_string(output_buffer, this->settings.custom_delimiter_char, (*objects.meta)[p], this->header);

			// Filter options
			if(settings.load_set_membership) this->printFILTER(output_buffer, p, objects);
			else output_buffer += '.';

			if(settings.load_info || settings.load_format || this->block.n_info_loaded) output_buffer += this->settings.custom_delimiter_char;
			else {
				output_buffer += '\n';
				continue;
			}

			// Print normal or with a custom delimiter
			this->printINFOVCF(output_buffer, p, objects);

			if(settings.load_format || this->block.n_format_loaded){
				output_buffer += this->settings.custom_delimiter_char;

				if(this->settings.format_ID_list.size())
					this->printFORMATCustom(output_buffer, this->settings.custom_delimiter_char, p, objects, genotypes_unpermuted);
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

		U32 info_match_limit = 1; // any match
		//info_match_limit = this->settings.info_list.size(); // all match

		U32 n_records_returned = 0;

		for(U32 p = 0; p < objects.meta->size(); ++p){
			//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
			//	continue;

			++n_records_returned;

			utility::to_vcf_string(output_buffer, this->settings.custom_delimiter_char, (*objects.meta)[p], this->header, this->settings.custom_output_controller);

			//output_buffer += '{';
			//utility::to_json(output_buffer, (*objects.meta)[p], this->header, this->settings.custom_output_controller);
			//output_buffer += '}';

			this->printINFOCustom(output_buffer, this->settings.custom_delimiter_char, p, objects);

			if(settings.load_format || this->block.n_format_loaded){
				output_buffer += this->settings.custom_delimiter_char;

				// Add FORMAT data
				this->printFORMATCustom(output_buffer, this->settings.custom_delimiter_char, p, objects, genotypes_unpermuted);
			}
			output_buffer += '\n';

			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

		return(n_records_returned);
	}

	/**<
	 *
	 * @param outputBuffer
	 * @param position
	 * @param objects
	 */
	void printFILTER(buffer_type& outputBuffer,
					 const U32& position,
					 const objects_type& objects) const
	{
		if(this->block.footer.n_filter_streams){
			const U32& n_filter_keys = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].n_keys;
			const U32* filter_keys   = this->block.footer.filter_bit_vectors[(*objects.meta)[position].filter_pattern_id].local_keys;
			if(n_filter_keys){
				// Local key -> global key
				outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
				for(U32 i = 1; i < n_filter_keys; ++i){
					outputBuffer += ';';
					outputBuffer += this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
				}
			} else
				outputBuffer += '.';
		} else {
			outputBuffer += '.';
		}
	}

	/**<
	 *
	 * @param buffer
	 * @param position
	 * @param objects
	 * @param genotypes_unpermuted
	 */
	void printFORMATVCF(buffer_type& buffer,
                        const U32& position,
						const objects_type& objects,
                        std::vector<core::GTObject>& genotypes_unpermuted) const
	{
		if(settings.load_format && this->block.n_format_loaded){
			if(this->block.n_format_loaded){
				const U32& n_format_keys = this->block.footer.format_bit_vectors[objects.meta->at(position).format_pattern_id].n_keys;
				const U32* format_keys   = this->block.footer.format_bit_vectors[objects.meta->at(position).format_pattern_id].local_keys;
				if(n_format_keys){
					// Print key map
					buffer += this->header.format_fields[this->block.footer.format_offsets[format_keys[0]].data_header.global_key].ID;
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						buffer += this->header.format_fields[this->block.footer.format_offsets[format_keys[i]].data_header.global_key].ID;
					}
					buffer += '\t';

					// Todo: print if no GT data
					// Begin print FORMAT data for each sample
					if(this->settings.load_ppa && this->block.header.controller.hasGTPermuted){
						objects.genotypes->at(position).getObjects(genotypes_unpermuted, this->header.getSampleNumber(), this->block.ppa_manager);
					} else {
						objects.genotypes->at(position).getObjects(genotypes_unpermuted, this->header.getSampleNumber());
					}

					buffer << genotypes_unpermuted[0];
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						objects.format_fields[format_keys[i]]->to_vcf_string(buffer, position, 0);
					}

					for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
						buffer += '\t';
						buffer << genotypes_unpermuted[s];
						for(U32 i = 1; i < n_format_keys; ++i){
							buffer  += ':';
							objects.format_fields[format_keys[i]]->to_vcf_string(buffer, position, s);
						}
					}
				} else // have no keys
					buffer += ".\t";
			}
		}
	}

	/**<
	 *
	 * @param outputBuffer
	 * @param delimiter
	 * @param position
	 * @param objects
	 * @param genotypes_unpermuted
	 */
	void printFORMATCustom(buffer_type& outputBuffer,
						   const char& delimiter,
						   const U32& position,
						   const objects_type& objects,
						   std::vector<core::GTObject>& genotypes_unpermuted) const
	{
		if(settings.load_format || this->block.n_format_loaded){
			const std::vector<U32>& targetKeys = objects.local_match_keychain_format[objects.meta->at(position).format_pattern_id];
			if(this->block.n_info_loaded && targetKeys.size()){
				// Print key map
				outputBuffer += this->header.format_fields[this->block.format_containers[targetKeys[0]].header.getGlobalKey()].ID;
				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer += ':';
					outputBuffer += this->header.format_fields[this->block.format_containers[targetKeys[i]].header.getGlobalKey()].ID;
				};

				outputBuffer += delimiter;

				// First individual
				//if(this->header.getSampleNumber() > 1) outputBuffer += '[';

				//if(targetKeys.size() > 1) outputBuffer += '[';
				objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, 0);
				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer += ':';
					objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, 0);
				}
				//if(targetKeys.size() > 1) outputBuffer += ']';

				for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
					outputBuffer += delimiter;
					//if(targetKeys.size() > 1) outputBuffer += '[';
					objects.format_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position, s);
					for(U32 i = 1; i < targetKeys.size(); ++i){
						outputBuffer  += ':';
						objects.format_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position, s);
					}
					//if(targetKeys.size() > 1) outputBuffer += ']';
				}


				//if(this->header.getSampleNumber() > 1) outputBuffer += ']';
			}
		}
	}

	/**<
	 * Outputs INFO fields in VCF formatting
	 * @param outputBuffer
	 * @param position
	 * @param objects
	 */
	void printINFOVCF(buffer_type& outputBuffer,
                      const U32& position,
					  const objects_type& objects) const
	{
		if(settings.load_info || this->block.n_info_loaded){
			const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta->at(position).info_pattern_id];

			if(this->block.n_info_loaded && targetKeys.size()){
				// First
				outputBuffer += objects.info_field_names[targetKeys[0]];
				if(objects.info_fields[targetKeys[0]]->emptyPosition(position) == false){
					outputBuffer += '=';
					objects.info_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position);
				}

				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer += ";";
					outputBuffer += objects.info_field_names[targetKeys[i]];
					if(this->header.info_fields[this->block.footer.info_offsets[targetKeys[i]].data_header.global_key].primitive_type == YON_VCF_HEADER_FLAG){
						continue;
					}
					if(objects.info_fields[targetKeys[i]]->emptyPosition(position)) continue;
					outputBuffer += '=';
					objects.info_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position);
				}
			} else
				outputBuffer += '.';
		}
	}

	/**<
	 * Outputs INFO fields in a custom format with the delimiter provided
	 * @param outputBuffer
	 * @param delimiter
	 * @param position
	 * @param objects
	 */
	void printINFOCustom(buffer_type& outputBuffer,
                         const char& delimiter,
						 const U32& position,
						 const objects_type& objects) const
	{
		if(settings.load_info || this->block.n_info_loaded){
			const std::vector<U32>& targetKeys = objects.local_match_keychain_info[objects.meta->at(position).info_pattern_id];
			if(this->block.n_info_loaded && targetKeys.size()){
				// Check if this target container is a FLAG
				if(this->header.info_fields[this->block.info_containers[targetKeys[0]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
					outputBuffer += this->header.info_fields[this->block.info_containers[targetKeys[0]].header.getGlobalKey()].ID;
				} else {
					// Check if the positon is empty
					if(objects.info_fields[targetKeys[0]]->emptyPosition(position) == false){
						objects.info_fields[targetKeys[0]]->to_vcf_string(outputBuffer, position);
					}
				}

				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer += delimiter;
					if(this->header.info_fields[this->block.info_containers[targetKeys[i]].header.getGlobalKey()].primitive_type == YON_VCF_HEADER_FLAG){
						outputBuffer += this->header.info_fields[this->block.info_containers[targetKeys[i]].header.getGlobalKey()].ID;
						continue;
					}
					if(objects.info_fields[targetKeys[i]]->emptyPosition(position)) continue;
					objects.info_fields[targetKeys[i]]->to_vcf_string(outputBuffer, position);
				}
			}
		}
	}


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
