#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "algorithm/timer.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/primitive_group_container.h"
#include "containers/meta_container.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "containers/variant_block.h"
#include "core/genotype_object.h"
#include "core/footer/footer.h"
#include "core/header/variant_header.h"
#include "math/fisher_math.h"
#include "math/square_matrix.h"
#include "math/basic_vector_math.h"
#include "utility/support_vcf.h"
#include "index/index.h"
#include "core/variant_reader_settings.h"
#include "core/variant_reader_objects.h"
#include "core/variant_reader_filters.h"

namespace tachyon{

class VariantReader{
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef core::VariantHeader                    header_type;
	typedef core::Footer                           footer_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef DataBlockSettings                      block_settings_type;
	typedef VariantReaderSettings                  settings_type;
	typedef index::Index                           index_type;
	typedef index::IndexEntry                      index_entry_type;
	typedef algorithm::VariantDigestManager checksum_type;
	typedef encryption::Keychain<>                 keychain_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;
	typedef containers::GenotypeSummary            genotype_summary_type;
	typedef VariantReaderFilters                   variant_filter_type;

	// Function pointers
	typedef void (self_type::*print_format_function)(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	typedef void (self_type::*print_info_function)(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	typedef void (self_type::*print_filter_function)(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	typedef buffer_type& (*print_meta_function)(buffer_type& buffer, const char& delimiter, const meta_entry_type& meta_entry, const header_type& header, const block_settings_type& controller);

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
	 * @return A reference instance of the block settings object
	 */
	inline block_settings_type& getBlockSettings(void){ return(this->block_settings); }

	/**<
	 * Retrieve current settings for the variant reader. This settings
	 * object controls the parsing/output of the reader itself. This is
	 * unlike the `block_settings_type` that controls the `DataBlock`
	 * parsing.
	 * @return A reference instance of the settings object
	 */
	inline settings_type& getSettings(void){ return(this->settings); }

	inline variant_filter_type& getFilterSettings(void){ return(this->variant_filters); }

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
		this->settings.input = filename;
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
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool getBlock(const index_entry_type& index_entry);

	/**<
	 * Get the current YON block in-order as a copy
	 * @return Returns a YON block. The container has a size of 0 upon fail/empty
	 */
	block_entry_type getBlock(void);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the block_settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool seek_to_block(const U32& blockID);

	/**<
	 *
	 * @return
	 */
	bool parseSettings(void);

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

		if(this->block_settings.annotate_extra){
			// fixme
			// if special
			// "FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"
			if(this->header.getInfoField("FS_A") == false)          this->header.literals += "\n##INFO=<ID=FS_A,Number=A,Type=Float>";
			if(this->header.getInfoField("AN") == false)            this->header.literals += "\n##INFO=<ID=AN,Number=A,Type=Integer>";
			if(this->header.getInfoField("NM") == false)            this->header.literals += "\n##INFO=<ID=NM,Number=A,Type=Integer>";
			if(this->header.getInfoField("NPM") == false)           this->header.literals += "\n##INFO=<ID=NPM,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC") == false)            this->header.literals += "\n##INFO=<ID=AC,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC_FWD") == false)        this->header.literals += "\n##INFO=<ID=AC_FWD,Number=A,Type=Integer>";
			if(this->header.getInfoField("AC_REV") == false)        this->header.literals += "\n##INFO=<ID=AC_REV,Number=A,Type=Integer>";
			if(this->header.getInfoField("HWE_P") == false)         this->header.literals += "\n##INFO=<ID=HWE_P,Number=A,Type=Float>";
			if(this->header.getInfoField("VT") == false)            this->header.literals += "\n##INFO=<ID=VT,Number=A,Type=String>";
			if(this->header.getInfoField("AF") == false)            this->header.literals += "\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">";
			if(this->header.getInfoField("MULTI_ALLELIC") == false) this->header.literals += "\n##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag>";
			if(this->header.getInfoField("F_PIC") == false)         this->header.literals += "\n##INFO=<ID=F_PIC,Number=A,Type=Float,Description=\"Population inbreeding coefficient (F-statistics)\">";
		}

		this->header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		this->header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
				  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

		this->header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + '\n';
		this->header.literals += this->getSettings().get_settings_string();

		// Output VCF header
		if(this->block_settings.show_vcf_header){
			this->header.writeVCFHeaderString(std::cout, this->block_settings.format_all.load || this->block_settings.format_list.size());
		}

		// If seek is active for targetted intervals
		if(this->getSettings().interval_strings.size()){
			if(this->getSettings().getBlockList().size()){
				for(U32 i = 0; i < this->getSettings().getBlockList().size(); ++i){
					if(this->getBlock(this->getSettings().getBlockList()[i]) == false){
						return(0);
					}
					n_variants += this->outputBlockVCF();
				}

				return(n_variants);
			} else { // Provided intervals but no matching YON blocks
				return(0);
			}
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
	const U32 outputBlockVCF(void){
		objects_type objects;
		this->loadObjects(objects);
		//this->variant_filters.filter_ref_allele("A", YON_CMP_EQUAL);
		//this->variant_filters.filter_alt_allele("[A]{2,}", YON_CMP_REGEX);
		this->variant_filters.build();

		// Reserve memory for output buffer
		// This is much faster than writing directly to ostream because of synchronisation
		io::BasicBuffer output_buffer(256000);
		if(this->block_settings.format_all.load) output_buffer.resize(256000 + this->header.getSampleNumber()*2);

		// Todo: in cases of non-diploid
		std::vector<core::GTObject> genotypes_unpermuted(this->header.getSampleNumber());
		for(U32 i = 0; i < this->header.getSampleNumber(); ++i){
			genotypes_unpermuted[i].alleles = new core::GTObjectAllele;
		}

		print_format_function print_format = &self_type::printFORMATDummy;
		if(this->block_settings.format_ID_list.size()) print_format = &self_type::printFORMATCustom;
		else if(block_settings.format_all.display)     print_format = &self_type::printFORMATVCF;
		print_info_function   print_info   = &self_type::printINFOVCF;
		print_meta_function   print_meta   = &utility::to_vcf_string;
		print_filter_function print_filter = &self_type::printFILTER;

		// Cycling over loaded meta objects
		for(U32 p = 0; p < objects.meta_container->size(); ++p){
			if(this->variant_filters.filter(objects, p) == false)
				continue;

			if(this->block_settings.custom_output_format)
				utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->header, this->block_settings);
			else
				utility::to_vcf_string(output_buffer, '\t', objects.meta_container->at(p), this->header);

			// Filter options
			(this->*print_filter)(output_buffer, p, objects);
			(this->*print_info)(output_buffer, '\t', p, objects);
			if(this->block_settings.annotate_extra) this->getGenotypeSummary(output_buffer, p, objects); // Todo: fixme
			(this->*print_format)(output_buffer, '\t', p, objects, genotypes_unpermuted);
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

		return(objects.meta_container->size());
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
		//U32 info_match_limit = 1; // any match
		//info_match_limit = this->block_settings.info_list.size(); // all match

		// Function pointer to use
		print_format_function print_format = &self_type::printFORMATDummy;
		print_info_function   print_info   = &self_type::printINFODummy;
		print_meta_function   print_meta   = &utility::to_vcf_string;
		print_filter_function print_filter = &self_type::printFILTERDummy;

		if(block_settings.output_json) print_meta = &utility::to_json_string;

		if(block_settings.format_all.display || this->block.n_format_loaded){
			if(block_settings.output_json){
				print_format = &self_type::printFORMATCustomVectorJSON;
			} else {
				if(block_settings.output_format_vector) print_format = &self_type::printFORMATCustomVector;
				else print_format = &self_type::printFORMATCustom;
			}
		}

		if(block_settings.info_all.display || this->block.n_info_loaded){
			if(block_settings.output_json) print_info = &self_type::printINFOCustomJSON;
			else print_info = &self_type::printINFOCustom;
		}

		if(block_settings.display_filter){
			if(block_settings.output_json) print_filter = &self_type::printFILTERJSON;
			else print_filter = &self_type::printFILTERCustom;
		}

		U32 n_records_returned = 0;

		if(block_settings.output_json) output_buffer += "\"block\":[";
		for(U32 position = 0; position < objects.meta_container->size(); ++position){
			//if(info_keep[objects.meta->at(p).getInfoPatternID()] < info_match_limit)
			//	continue;

			if(block_settings.output_json){
				if(position != 0) output_buffer += ",\n";
				output_buffer += "{";
			}
			++n_records_returned;

			(*print_meta)(output_buffer, this->block_settings.custom_delimiter_char, objects.meta_container->at(position), this->header, this->block_settings);
			(this->*print_filter)(output_buffer, position, objects);
			(this->*print_info)(output_buffer, this->block_settings.custom_delimiter_char, position, objects);
			(this->*print_format)(output_buffer, this->block_settings.custom_delimiter_char, position, objects, genotypes_unpermuted);

			if(block_settings.output_json) output_buffer += "}";
			else output_buffer += '\n';
			//output_buffer += "}";

			// Flush if buffer is large
			if(output_buffer.size() > 65536){
				std::cout.write(output_buffer.data(), output_buffer.size());
				output_buffer.reset();
				std::cout.flush();
			}
		}
		if(block_settings.output_json) output_buffer += "]";

		// Flush buffer
		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		std::cout.flush();

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

	// Calculations
	TACHYON_VARIANT_CLASSIFICATION_TYPE classifyVariant(const meta_entry_type& meta, const U32& allele) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool parseIntervals(void){
		// Intervals pass expression tests
		if(this->getSettings().validateIntervalStrings() == false)
			return(false);

		// Data reference
		std::vector<std::string>& intervals = this->getSettings().interval_strings;

		// No intervals to parse
		if(intervals.size() == 0)
			return true;

		// Blocks to load
		std::vector<index_entry_type>& blocks_to_load = this->getSettings().getBlockList();

		// Parse each interval
		for(U32 i = 0; i < intervals.size(); ++i){
			intervals[i] = utility::remove_whitespace(intervals[i]);
			core::HeaderContig* contig = nullptr;

			if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_ONLY )){
				std::cerr << "chromosome only" << std::endl;
				if(!this->header.getContig(intervals[i],contig)){
					std::cerr << "cant find contig: " << intervals[i] << std::endl;
					return(false);
				}

				std::cerr << "Parsed: " << intervals[i] << std::endl;

			} else if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_POSITION )){
				std::cerr << "chromosome pos" << std::endl;
				std::vector<std::string> substrings = utility::split(intervals[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!this->header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				U64 position = atof(substrings[1].data());
				std::cerr << "Parsed: " << substrings[0] << "," << position << std::endl;

				std::vector<index_entry_type> target_blocks = this->index.findOverlap(contig->contigID, position);
				blocks_to_load.insert( blocks_to_load.end(), target_blocks.begin(), target_blocks.end() );

			} else if (std::regex_match (intervals[i], constants::YON_REGEX_CONTIG_RANGE )){
				std::cerr << "chromosome pos - pos" << std::endl;
				std::vector<std::string> substrings = utility::split(intervals[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!this->header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				std::vector<std::string> position_strings = utility::split(substrings[1], '-');
				if(position_strings[0].size() == 0 || position_strings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}
				U64 position_from = atof(position_strings[0].data());
				U64 position_to   = atof(position_strings[1].data());
				if(position_from > position_to) std::swap(position_from, position_to);

				std::cerr << "Parsed: " << substrings[0] << "," << position_from << "," << position_to << std::endl;

				std::vector<index_entry_type> target_blocks = this->index.findOverlap(contig->contigID, position_from, position_to);
				blocks_to_load.insert( blocks_to_load.end(), target_blocks.begin(), target_blocks.end() );

			} else {
				std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << intervals[i] << std::endl;
				return false;
			}
		}

		if(blocks_to_load.size() == 0)
			return true;

		// Dedupe
		std::sort(blocks_to_load.begin(), blocks_to_load.end());
		std::vector<index_entry_type> blocks_to_load_deduped;
		blocks_to_load_deduped.push_back(blocks_to_load[0]);
		for(U32 i = 1; i < blocks_to_load.size(); ++i){
			if(blocks_to_load_deduped.back() != blocks_to_load[i]){
				blocks_to_load_deduped.push_back(blocks_to_load[i]);
			}
		}

		// Debug
		//std::cerr << "size: " << blocks_to_load_deduped.size() << std::endl;
		//for(U32 i = 0; i < blocks_to_load_deduped.size(); ++i){
		//	std::cerr << i << "\t";
		//	blocks_to_load_deduped[i].print(std::cerr);
		//	std::cerr << std::endl;
		//}

		// Assign deduped vector to settings reference
		blocks_to_load = blocks_to_load_deduped;

		return true;
	}


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 timings_meta(){
		containers::MetaContainer meta(this->block);
		buffer_type temp(meta.size() * 1000);

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(temp, this->block_settings.custom_delimiter_char, meta[p], this->header);
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

			std::cout << (int)objects_true[i].alleles[0].allele << (objects_true[i].alleles[1].phase ? '/' : '|') << (int)objects_true[i].alleles[1].allele;
			for(U32 i = 1; i < objects_true.size(); ++i){
				std::cout << '\t' << (int)objects_true[i].alleles[0].allele << (objects_true[i].alleles[1].phase ? '/' : '|') << (int)objects_true[i].alleles[1].allele;
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

	std::vector<double> calculateStrandBiasAlleles(const meta_entry_type& meta, const genotype_summary_type& genotype_summary, const bool phred_scale = true) const{
		std::vector<double> strand_bias_p_values(meta.n_alleles);
		double fisher_left_p, fisher_right_p, fisher_twosided_p;

		kt_fisher_exact(
		genotype_summary.vectorA_[2], // A: Allele on forward strand
		genotype_summary.vectorB_[2], // B: Allele on reverse strand
		genotype_summary.alleleCountA() - (genotype_summary.vectorA_[2]), // C: Not allele on forward strand
		genotype_summary.alleleCountB() - (genotype_summary.vectorB_[2]), // D: Not allele on reverse strand
		&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

		if(phred_scale) strand_bias_p_values[0] = std::abs(-10 * log10(fisher_twosided_p));
		else strand_bias_p_values[0] = fisher_twosided_p;

		// If n_alleles = 2 then they are identical because of symmetry
		if(meta.n_alleles > 2){
			for(U32 p = 1; p < meta.n_alleles; ++p){
				kt_fisher_exact(
				genotype_summary.vectorA_[2+p], // A: Allele on forward strand
				genotype_summary.vectorB_[2+p], // B: Allele on reverse strand
				genotype_summary.alleleCountA() - (genotype_summary.vectorA_[2+p]), // C: Not allele on forward strand
				genotype_summary.alleleCountB() - (genotype_summary.vectorB_[2+p]), // D: Not allele on reverse strand
				&fisher_left_p, &fisher_right_p, &fisher_twosided_p);

				if(phred_scale) strand_bias_p_values[p] = std::abs(-10 * log10(fisher_twosided_p));
				else strand_bias_p_values[p] = fisher_twosided_p;
			}
		}
		return(strand_bias_p_values);
	}

	void getGenotypeSummary(buffer_type& buffer, const U32& position, objects_type& objects) const{
		if(this->block_settings.alleles.load == false || this->block_settings.genotypes_all.load == false || this->block_settings.controller.load == false || this->block_settings.set_membership.load == false){
			std::cerr << utility::timestamp("ERROR") << "Cannot run function without loading: SET-MEMBERSHIP, GT, REF or ALT, CONTIG or POSITION..." << std::endl;
			return;
		}

		if(buffer.back() != ';' && buffer.back() != '\t') buffer += ';';

		//objects_type objects;
		//this->loadObjects(objects);
		//U32 n_variants_parsed = 0;

		//for(U32 i = 0; i < objects.genotypes->size(); ++i){
			if(objects.meta_container->at(position).isDiploid() == false){
				std::cerr << "is not diploid" << std::endl;
				return;
			}

			// If set membership is -1 then calculate all fields
			// Set target FLAG set to all ones; update with actual values if they exist
			U16 target_flag_set = 65535;
			if(objects.meta_container->at(position).getInfoPatternID() != -1)
				target_flag_set = objects.additional_info_execute_flag_set[objects.meta_container->at(position).getInfoPatternID()];

			// Get genotype summary data
			objects.genotype_container->at(position).getSummary(*objects.genotype_summary);
			std::vector<double> hwe_p = objects.genotype_summary->calculateHardyWeinberg(objects.meta_container->at(position));
			std::vector<double> af    = objects.genotype_summary->calculateAlleleFrequency(objects.meta_container->at(position));


			//utility::to_vcf_string(stream, this->block_settings.custom_delimiter_char, meta, this->header);

			if(target_flag_set & 1){
				std::vector<double> allele_bias = this->calculateStrandBiasAlleles(objects.meta_container->at(position), *objects.genotype_summary, true);
				buffer += "FS_A=";
				buffer.AddReadble(allele_bias[0]);
				for(U32 p = 1; p < allele_bias.size(); ++p){
					buffer += ',';
					buffer.AddReadble(allele_bias[p]);
				}
			}

			if(target_flag_set & 2){
				buffer += ";AN=";
				buffer.AddReadble(objects.genotype_summary->alleleCount());
			}

			if(target_flag_set & 4){
				if(objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1]){
					buffer += ";NM=";
					buffer.AddReadble(objects.genotype_summary->vectorA_[1] + objects.genotype_summary->vectorB_[1]);
				}
			}

			if(target_flag_set & 8){
				if(objects.genotype_summary->vectorA_[0] + objects.genotype_summary->vectorB_[0]){
					buffer += ";NPM=";
					buffer.AddReadble(objects.genotype_summary->vectorA_[0] + objects.genotype_summary->vectorB_[0]);
				}
			}

			if(target_flag_set & 16){
				buffer += ";AC=";
				buffer.AddReadble(objects.genotype_summary->vectorA_[2] + objects.genotype_summary->vectorB_[2]);
				for(U32 p = 1; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorA_[2+p] + objects.genotype_summary->vectorB_[2+p]);
				}
			}

			if(target_flag_set & 32){
				buffer += ";AC_FWD=";
				buffer.AddReadble(objects.genotype_summary->vectorA_[2]);
				for(U32 p = 1; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorA_[2+p]);
				}
			}

			if(target_flag_set & 64){
				buffer += ";AC_REV=";
				buffer.AddReadble(objects.genotype_summary->vectorB_[2]);
				for(U32 p = 1; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ",";
					buffer.AddReadble(objects.genotype_summary->vectorB_[2+p]);
				}
			}

			if(target_flag_set & 128){
				buffer += ";AF=";
				buffer.AddReadble(af[0]);
				for(U32 p = 1; p < af.size(); ++p){
					buffer += ",";
					buffer.AddReadble(af[p]);
				}
			}

			if(target_flag_set & 256){
				buffer += ";HWE_P=";
				buffer.AddReadble(hwe_p[0]);
				for(U32 p = 1; p < hwe_p.size(); ++p){
					buffer += ",";
					buffer.AddReadble(hwe_p[p]);
				}
			}

			if(target_flag_set & 512){
				// Classify
				buffer += ";VT=";
				buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->classifyVariant(objects.meta_container->at(position), 1)];

				for(U32 p = 2; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ',';
					buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->classifyVariant(objects.meta_container->at(position), p)];
				}
			}

			if(target_flag_set & 1024){
				if(objects.meta_container->at(position).n_alleles > 2) buffer += ";MULTI_ALLELIC";
			}

			// Population inbreeding coefficient: F = (Hexp - Hobs) /Hexp
			if(target_flag_set & 2048){
				// Allele frequency of A
				const double p = ((double)2*objects.genotype_summary->matrix_[2][2] + objects.genotype_summary->matrix_[2][3] + objects.genotype_summary->matrix_[3][2]) / (2*objects.genotype_summary->genotypeCount());
				// Genotype frequency of heterozyotes
				const double pg = ((double)objects.genotype_summary->matrix_[2][3] + objects.genotype_summary->matrix_[3][2]) / objects.genotype_summary->genotypeCount();
				// Expected heterozygosity
				const double exp = 2*p*(1-p);
				// Population inbreeding coefficient: F
				const double f_pic = exp > 0 ? (exp-pg)/exp : 0;
				buffer += ";F_PIC=";
				buffer.AddReadble(f_pic);
			}

			objects.genotype_summary->clear();
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
	std::ifstream       stream;
	U64                 filesize;

	// Actual data
	block_entry_type    block;

	// Supportive objects
	block_settings_type block_settings;
	settings_type       settings;
	variant_filter_type variant_filters;

	header_type         header;
	footer_type         footer;
	index_type          index;
	checksum_type       checksums;
	codec_manager_type  codec_manager;
	keychain_type       keychain;
};

}

#endif /* CORE_TACHYON_READER_H_ */
