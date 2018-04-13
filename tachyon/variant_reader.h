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

class VariantReader{
	typedef VariantReader                          self_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef core::VariantHeader                    header_type;
	typedef core::Footer                           footer_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef core::DataBlockSettings                settings_type;
	typedef index::Index                           index_type;
	typedef algorithm::VariantDigitalDigestManager checksum_type;
	typedef encryption::Keychain                   keychain_type;
	typedef core::MetaEntry                        meta_entry_type;

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

	bool parseSettings(void){
		settings.load_info_ID_loaded.clear();
		settings.load_format_ID_loaded.clear();

		// Map INFO
		U32 info_matches = 0;
		for(U32 i = 0; i < settings.info_list.size(); ++i){
			const S32 global_key = this->has_info_field(settings.info_list[i]);
			if(global_key){
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
						//std::cerr << "match @ local: " << i << std::endl;
						break;
					}
				}

				if(local_key == -1){
					//std::cerr << "could not find local" << std::endl;
				}
			}
		}

		// Todo: for each local hit
		// build new bit-vector

		return(true);
	}


	/**<
	 *
	 * @param stream
	 * @return
	 */
	U64 outputVCF(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer* gt = nullptr;
		if(this->block.header.controller.hasGT && settings.loadGenotypesRLE_)
			gt = new containers::GenotypeContainer(this->block, meta);

		containers::FormatContainerInterface** fits = new containers::FormatContainerInterface*[this->block.n_format_loaded];

		if(this->block.n_format_loaded){
			for(U32 i = 0; i < this->block.n_format_loaded; ++i){
				//std::cerr << this->block.format_containers[i].buffer_data_uncompressed.size() << std::endl;
				const U32 global_key = settings.load_format_ID_loaded[i].offset->data_header.global_key;
				//std::cerr << i << "/" << this->block.footer.n_format_streams << ": " << global_key << "@" << this->header.format_fields[this->block.footer.format_offsets[i].data_header.global_key].ID << '\t' << this->header.format_fields[global_key].getType() << std::endl;
				std::vector<bool> matches = this->get_format_field_pattern_matches(this->header.format_fields[global_key].ID);

				if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					fits[i] = new containers::FormatContainer<S32>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());
					//global_fields.push_back(this->header.format_fields[global_key].ID);
					//std::cerr << "s32: " << fits[i]->size() << std::endl;
				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					fits[i] = new containers::FormatContainer<std::string>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());
					//global_fields.push_back(this->header.format_fields[global_key].ID);
					//std::cerr << "string: " << fits[i]->size() << std::endl;

				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					fits[i] = new containers::FormatContainer<float>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());
					//global_fields.push_back(this->header.format_fields[global_key].ID);
					//std::cerr << "float: " << fits[i]->size() << std::endl;
				} else {
					fits[i] = new containers::FormatContainer<U32>;
					//global_fields.push_back(this->header.format_fields[global_key].ID);
				}
			}
		}

		// Store as double pointers to avoid memory collisions because
		// info containers have different class members
		containers::InfoContainerInterface** its = new containers::InfoContainerInterface*[this->block.n_info_loaded];

		std::vector<std::string> global_fields;
		if(this->block.n_info_loaded){
			for(U32 i = 0; i < this->block.n_info_loaded; ++i){
				const U32 global_key = settings.load_info_ID_loaded[i].offset->data_header.global_key;

				//std::cerr << i << "/" << this->block.n_info_loaded << ": " << global_key << "@" << this->header.info_fields[this->block.footer.info_offsets[i].data_header.global_key].ID << std::endl;
				std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[global_key].ID);

				if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					its[i] = new containers::InfoContainer<S32>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					its[i] = new containers::InfoContainer<std::string>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					its[i] = new containers::InfoContainer<float>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else {
					its[i] = new containers::InfoContainer<U32>();
					global_fields.push_back(this->header.info_fields[global_key].ID);
				}
			}
		}

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(stream, this->settings.customDelimiterChar_, meta[p], this->header);

			if(settings.loadSetMembership_){
				if(this->block.footer.n_filter_streams){
					const U32& n_filter_keys = this->block.footer.filter_bit_vectors[meta[p].filter_pattern_id].n_keys;
					const U32* filter_keys   = this->block.footer.filter_bit_vectors[meta[p].filter_pattern_id].local_keys;
					if(n_filter_keys){
						// Local key -> global key
						stream << this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[0]].data_header.global_key].ID;
						for(U32 i = 1; i < n_filter_keys; ++i){
							stream << ';' << this->header.filter_fields[this->block.footer.filter_offsets[filter_keys[i]].data_header.global_key].ID;
						}
					} else
						stream.put('.');
				} else {
					stream.put('.');
				}
			} else
				stream.put('.');

			if(settings.loadINFO_ || settings.loadFORMAT_ || this->block.n_info_loaded) stream.put('\t');
			else {
				stream.put('\n');
				continue;
			}

			// At this point we need to ascertain that the target set membership identifier
			// contain a loaded INFO value
			if(settings.loadINFO_ || this->block.n_info_loaded){
				if(this->block.n_info_loaded){
					const U32& n_keys = this->block.footer.info_bit_vectors[meta[p].info_pattern_id].n_keys;
					const U32* keys   = this->block.footer.info_bit_vectors[meta[p].info_pattern_id].local_keys;

					//for(U32 z = 0; z < this->settings.load_info_ID_loaded.size(); ++z){
					//	std::cerr << "Key: " << this->settings.load_info_ID_loaded[z].offset->data_header.global_key << " of " << this->settings.load_info_ID_loaded.size() << " " << this->block.footer.info_bit_vectors[meta[p].info_pattern_id][settings.load_info_ID_loaded[z].offset->data_header.global_key] << std::endl;
					//}

					// First
					stream << global_fields[keys[0]];
					if(its[keys[0]]->emptyPosition(p) == false){;
						stream.put('=');
						its[keys[0]]->to_vcf_string(stream, p);
					}

					for(U32 i = 1; i < n_keys; ++i){
						stream << ";" << global_fields[keys[i]];
						if(this->header.info_fields[this->block.footer.info_offsets[keys[i]].data_header.global_key].primitive_type == 2){
							continue;
						}
						if(its[keys[i]]->emptyPosition(p)) continue;
						stream.put('=');
						its[keys[i]]->to_vcf_string(stream, p);
					}
				} else
					stream.put('.');
			}

			if(settings.loadFORMAT_ && this->block.n_format_loaded) stream.put('\t');
			else stream.put('\n');

			this->printFORMAT(stream, meta[p], p, fits, gt);
		}

		if(settings.loadINFO_){
			for(U32 i = 0; i < this->block.n_info_loaded; ++i)
				delete its[i];
		}

		if(settings.loadFORMAT_){
			for(U32 i = 0; i < this->block.n_format_loaded; ++i)
				delete fits[i];
		}

		delete [] its; delete [] fits;
		delete gt;
		return(meta.size());
	}

	/**<
	 * This function assumes that ALL available FORMAT fields described in the header
	 * has been decrypted, decompressed, and loaded and are available for use
	 * @param stream
	 * @param meta_entry
	 * @param individual
	 * @param fits
	 * @param gt
	 */
	void printFORMAT(std::ostream& stream, const meta_entry_type& meta_entry, const U32& individual, containers::FormatContainerInterface** fits, containers::GenotypeContainer* gt){
		if(settings.loadFORMAT_ && this->block.n_format_loaded){
			if(this->block.n_format_loaded){
				std::vector<core::GTObject> objects_true(this->header.getSampleNumber());

				const U32& n_format_keys = this->block.footer.format_bit_vectors[meta_entry.format_pattern_id].n_keys;
				const U32* format_keys   = this->block.footer.format_bit_vectors[meta_entry.format_pattern_id].local_keys;
				if(n_format_keys){
					// Print key map
					stream << this->header.format_fields[this->block.footer.format_offsets[format_keys[0]].data_header.global_key].ID;
					for(U32 i = 1; i < n_format_keys; ++i){
						stream.put(':');
						stream << this->header.format_fields[this->block.footer.format_offsets[format_keys[i]].data_header.global_key].ID;
					}
					stream.put('\t');

					// Todo: print if no GT data
					// Begin print FORMAT data for each sample
					if(this->settings.loadPPA_ && this->block.header.controller.hasGTPermuted){
						gt->at(individual).getObjects(objects_true, this->header.getSampleNumber(), this->block.ppa_manager);
					} else {
						gt->at(individual).getObjects(objects_true, this->header.getSampleNumber());
					}

					stream << objects_true[0];
					for(U32 i = 1; i < n_format_keys; ++i){
						stream.put(':');
						fits[format_keys[i]]->to_vcf_string(stream, individual, 0);
					}

					for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
						stream.put('\t');
						stream << objects_true[s];
						for(U32 i = 1; i < n_format_keys; ++i){
							stream.put(':');
							fits[format_keys[i]]->to_vcf_string(stream, individual, s);
						}
					}
					stream.put('\n');
				} else // have no keys
					stream << ".\t";
			}
		}
	}

	U64 outputVCFBuffer(){
		containers::MetaContainer meta(this->block);
		containers::GenotypeContainer* gt = nullptr;
		if(this->block.header.controller.hasGT && settings.loadGenotypesRLE_)
			gt = new containers::GenotypeContainer(this->block, meta);

		containers::FormatContainerInterface** fits = new containers::FormatContainerInterface*[this->block.n_format_loaded];

		if(this->block.n_format_loaded){
			for(U32 i = 0; i < this->block.n_format_loaded; ++i){
				const U32 global_key = settings.load_format_ID_loaded[i].offset->data_header.global_key;
				std::vector<bool> matches = this->get_format_field_pattern_matches(this->header.format_fields[global_key].ID);

				if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					fits[i] = new containers::FormatContainer<S32>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());
				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.format_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					fits[i] = new containers::FormatContainer<std::string>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());

				} else if(this->header.format_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					fits[i] = new containers::FormatContainer<float>(this->block.format_containers[i], meta, matches, this->header.getSampleNumber());
				} else {
					fits[i] = new containers::FormatContainer<U32>;
				}
			}
		}

		// Store as double pointers to avoid memory collisions because
		// info containers have different class members
		containers::InfoContainerInterface** its = new containers::InfoContainerInterface*[this->block.n_info_loaded];

		std::vector<std::string> global_fields;
		if(this->block.n_info_loaded){
			for(U32 i = 0; i < this->block.n_info_loaded; ++i){
				const U32 global_key = settings.load_info_ID_loaded[i].offset->data_header.global_key;
				std::vector<bool> matches = this->get_info_field_pattern_matches(this->header.info_fields[global_key].ID);

				if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_INTEGER){
					its[i] = new containers::InfoContainer<S32>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_STRING ||
						  this->header.info_fields[global_key].getType() == YON_VCF_HEADER_CHARACTER){
					its[i] = new containers::InfoContainer<std::string>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else if(this->header.info_fields[global_key].getType() == YON_VCF_HEADER_FLOAT){
					its[i] = new containers::InfoContainer<float>(this->block.info_containers[i], meta, matches);
					global_fields.push_back(this->header.info_fields[global_key].ID);
				} else {
					its[i] = new containers::InfoContainer<U32>();
					global_fields.push_back(this->header.info_fields[global_key].ID);
				}
			}
		}

		// Todo: Efficiency please!
		std::vector< std::vector<U32> > local_match_keychain_info(this->block.footer.n_info_patterns);

		// If loading all INFO values then return them in the ORIGINAL order
		if(this->settings.loadINFO_){
			for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // Number of info patterns
				for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
					for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
						if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
							local_match_keychain_info[i].push_back(k);
						}
					}
				}
			}
		}
		// If loading custom INFO fields then return them in the REQUESTED order
		else {
			for(U32 i = 0; i < this->block.footer.n_info_patterns; ++i){ // Number of info patterns
				for(U32 k = 0; k < settings.load_info_ID_loaded.size(); ++k){ // Number of loaded INFO identifiers
					for(U32 j = 0; j < this->block.footer.info_bit_vectors[i].n_keys; ++j){ // Number of keys in pattern [i]
						if(this->block.footer.info_offsets[this->block.footer.info_bit_vectors[i].local_keys[j]].data_header.global_key == settings.load_info_ID_loaded[k].offset->data_header.global_key){
							local_match_keychain_info[i].push_back(k);
						}
					}
				}
			}
		}

		io::BasicBuffer outputBuffer(256000 + this->header.getSampleNumber()*2);
		std::vector<core::GTObject> objects_true(this->header.getSampleNumber());

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(outputBuffer, this->settings.customDelimiterChar_, meta[p], this->header);

			if(settings.loadSetMembership_){
				if(this->block.footer.n_filter_streams){
					const U32& n_filter_keys = this->block.footer.filter_bit_vectors[meta[p].filter_pattern_id].n_keys;
					const U32* filter_keys   = this->block.footer.filter_bit_vectors[meta[p].filter_pattern_id].local_keys;
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
			} else
				outputBuffer += '.';

			if(settings.loadINFO_ || settings.loadFORMAT_ || this->block.n_info_loaded) outputBuffer += this->settings.customDelimiterChar_;
			else {
				outputBuffer += '\n';
				continue;
			}

			// Print normal or with a custom delimiter
			if(this->settings.customDelimiter_)
				this->printINFODelimited(outputBuffer, this->settings.customDelimiterChar_, meta[p], p, its, global_fields, local_match_keychain_info);
			else
				this->printINFO(outputBuffer, meta[p], p, its, global_fields, local_match_keychain_info);

			if(settings.loadFORMAT_ && this->block.n_format_loaded) outputBuffer += this->settings.customDelimiterChar_;
			else outputBuffer += '\n';

			// Add FORMAT data
			this->printFORMAT(outputBuffer, meta[p], p, fits, gt, objects_true);

			if(outputBuffer.size() > 65536){
				std::cout.write(outputBuffer.data(), outputBuffer.size());
				outputBuffer.reset();
				std::cout.flush();
			}
		}

		std::cout.write(outputBuffer.data(), outputBuffer.size());
		outputBuffer.reset();
		std::cout.flush();

		if(settings.loadINFO_){
			for(U32 i = 0; i < this->block.n_info_loaded; ++i)
				delete its[i];
		}

		if(settings.loadFORMAT_){
			for(U32 i = 0; i < this->block.n_format_loaded; ++i)
				delete fits[i];
		}

		delete [] its; delete [] fits;
		delete gt;
		return(meta.size());
	}

	void printFORMAT(buffer_type& buffer, const meta_entry_type& meta_entry, const U32& individual, containers::FormatContainerInterface** fits, containers::GenotypeContainer* gt, std::vector<core::GTObject>& objects_true) const{
		if(settings.loadFORMAT_ && this->block.n_format_loaded){
			if(this->block.n_format_loaded){

				const U32& n_format_keys = this->block.footer.format_bit_vectors[meta_entry.format_pattern_id].n_keys;
				const U32* format_keys   = this->block.footer.format_bit_vectors[meta_entry.format_pattern_id].local_keys;
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
					if(this->settings.loadPPA_ && this->block.header.controller.hasGTPermuted){
						gt->at(individual).getObjects(objects_true, this->header.getSampleNumber(), this->block.ppa_manager);
					} else {
						gt->at(individual).getObjects(objects_true, this->header.getSampleNumber());
					}

					buffer << objects_true[0];
					for(U32 i = 1; i < n_format_keys; ++i){
						buffer += ':';
						fits[format_keys[i]]->to_vcf_string(buffer, individual, 0);
					}

					for(U64 s = 1; s < this->header.getSampleNumber(); ++s){
						buffer += '\t';
						buffer << objects_true[s];
						for(U32 i = 1; i < n_format_keys; ++i){
							buffer  += ':';
							fits[format_keys[i]]->to_vcf_string(buffer, individual, s);
						}
					}
					buffer += '\n';
				} else // have no keys
					buffer += ".\t";
			}
		}
	}

	void printINFO(buffer_type& outputBuffer, const meta_entry_type& meta_entry, const U32& site, containers::InfoContainerInterface** its, const std::vector<std::string>& global_fields, const std::vector< std::vector<U32> >& local_match_keychain_info) const{
		if(settings.loadINFO_ || this->block.n_info_loaded){
			const std::vector<U32>& targetKeys = local_match_keychain_info[meta_entry.info_pattern_id];
			if(this->block.n_info_loaded && targetKeys.size()){

				// First
				outputBuffer += global_fields[targetKeys[0]];
				if(its[targetKeys[0]]->emptyPosition(site) == false){
					outputBuffer += '=';
					its[targetKeys[0]]->to_vcf_string(outputBuffer, site);
				}

				for(U32 i = 1; i < targetKeys.size(); ++i){
					outputBuffer += ";";
					outputBuffer += global_fields[targetKeys[i]];
					if(this->header.info_fields[this->block.footer.info_offsets[targetKeys[i]].data_header.global_key].primitive_type == YON_VCF_HEADER_FLAG){
						continue;
					}
					if(its[targetKeys[i]]->emptyPosition(site)) continue;
					outputBuffer += '=';
					its[targetKeys[i]]->to_vcf_string(outputBuffer, site);
				}
			} else
				outputBuffer += '.';
		}
	}

	void printINFODelimited(buffer_type& outputBuffer, const char& delimiter, const meta_entry_type& meta_entry, const U32& site, containers::InfoContainerInterface** its, const std::vector<std::string>& global_fields, const std::vector< std::vector<U32> >& local_match_keychain_info) const{
			if(settings.loadINFO_ || this->block.n_info_loaded){
				const std::vector<U32>& targetKeys = local_match_keychain_info[meta_entry.info_pattern_id];
				if(this->block.n_info_loaded && targetKeys.size()){

					// First
					//outputBuffer += global_fields[targetKeys[0]];
					if(its[targetKeys[0]]->emptyPosition(site) == false){
						//outputBuffer += '=';
						its[targetKeys[0]]->to_vcf_string(outputBuffer, site);
					}

					for(U32 i = 1; i < targetKeys.size(); ++i){
						//outputBuffer += ";";
						outputBuffer += delimiter;
						//outputBuffer += global_fields[targetKeys[i]];
						if(this->header.info_fields[this->block.footer.info_offsets[targetKeys[i]].data_header.global_key].primitive_type == YON_VCF_HEADER_FLAG){
							continue;
						}
						if(its[targetKeys[i]]->emptyPosition(site)) continue;
						//outputBuffer += '=';
						its[targetKeys[i]]->to_vcf_string(outputBuffer, site);
					}
				} else {
					//outputBuffer += '.';
				}
			}
		}


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 timings_meta(){
		containers::MetaContainer meta(this->block);
		buffer_type temp(meta.size() * 1000);

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(temp, this->settings.customDelimiterChar_, meta[p], this->header);
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
		//std::cerr << "before meta: " << this->block.size() << std::endl;
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

			//std::cerr << objects.size() << ' ' << objects_all.size() << std::endl;
			//for(U32 i = 0; i < objects.size(); ++i){
			//	std::cerr << (int)objects[i].alleles[0].first << "|" << (int)objects[i].alleles[1].first << "," << objects[i].n_objects << ' ';
			//}
			//std::cerr << std::endl;


			std::cout << (int)objects_true[i].alleles[0].first << (objects_true[i].alleles[1].second ? '/' : '|') << (int)objects_true[i].alleles[1].first;
			for(U32 i = 1; i < objects_true.size(); ++i){
				std::cout << '\t' << (int)objects_true[i].alleles[0].first << (objects_true[i].alleles[1].second ? '/' : '|') << (int)objects_true[i].alleles[1].first;
			}
			std::cout << std::endl;


			// Print the difference
			//std::cerr << objects.size() << '\t' << objects_all.size() << '\t' << objects_true.size() << std::endl;
			// Dump data
			//gt[i].getMeta().toVCFString(std::cout, this->header, this->block.header.contigID, this->block.header.minPosition);
			//utility::to_vcf_string(std::cout, objects_true) << '\n';
		}

		//containers::MetaContainer meta(this->block);
		//for(U32 i = 0; i < this->block.index_entry.n_format_streams; ++i)
		//	std::cerr << this->block.format_containers[i].buffer_data_uncompressed.size() << std::endl;
		//containers::InfoContainer<std::string>* it2 = this->get_balanced_info_container<std::string>("CSQ",meta);
		/*
		containers::InfoContainer<U32>* it2 = this->get_balanced_info_container<U32>("DP", meta);
		if(it2 != nullptr){
			std::cerr << it2->size() << std::endl;
		}
		delete it2;
		*/
		return(gt.size());
		//std::cerr << meta.size() << std::endl;
		//return(meta.size());
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
				utility::to_vcf_string(stream, this->settings.customDelimiterChar_, meta[i], this->header);
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
