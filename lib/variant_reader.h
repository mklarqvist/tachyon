#ifndef CORE_TACHYON_READER_H_
#define CORE_TACHYON_READER_H_

#include <regex>
#include <cmath>

#include "zstd.h"
#include "zstd_errors.h"

#include "algorithm/compression/compression_manager.h"
#include "algorithm/digest/variant_digest_manager.h"
#include "algorithm/encryption/encryption_decorator.h"
#include "algorithm/timer.h"
#include "containers/format_container.h"
#include "containers/format_container_string.h"
#include "containers/genotype_container.h"
#include "containers/info_container.h"
#include "containers/info_container_string.h"
#include "containers/interval_container.h"
#include "containers/meta_container.h"
#include "containers/primitive_group_container.h"
#include "containers/variant_block.h"
#include "containers/variant_block_container.h"
#include "core/footer/footer.h"
#include "core/genotype_object.h"
#include "core/variant_site_annotation.h"
#include "core/header/variant_header.h"
#include "core/variant_reader_filters.h"
#include "core/variant_reader_objects.h"
#include "core/variant_reader_settings.h"
#include "core/data_block_settings.h"
#include "index/index.h"
#include "math/basic_vector_math.h"
#include "math/fisher_math.h"
#include "math/square_matrix.h"
#include "utility/support_vcf.h"
#include "io/basic_reader.h"

namespace tachyon{

class VariantReader{
private:
	typedef VariantReader                          self_type;
	typedef io::BasicBuffer                        buffer_type;
	typedef VariantHeader                          header_type;
	typedef core::Footer                           footer_type;
	typedef core::MetaEntry                        meta_entry_type;
	typedef algorithm::CompressionManager          codec_manager_type;
	typedef DataBlockSettings                      block_settings_type;
	typedef VariantReaderSettings                  settings_type;
	typedef index::Index                           index_type;
	typedef index::IndexEntry                      index_entry_type;
	typedef algorithm::VariantDigestManager        checksum_type;
	typedef encryption::Keychain<>                 keychain_type;
	typedef VariantReaderObjects                   objects_type;
	typedef containers::VariantBlock               block_entry_type;
	typedef containers::MetaContainer              meta_container_type;
	typedef containers::GenotypeContainer          gt_container_type;
	typedef containers::InfoContainerInterface     info_interface_type;
	typedef containers::FormatContainerInterface   format_interface_type;
	typedef containers::IntervalContainer          interval_container_type;
	typedef containers::VariantBlockContainer      variant_container_type;
	typedef containers::GenotypeSummary            genotype_summary_type;
	typedef containers::VariantSiteAnnotation      site_annotation_type;
	typedef VariantReaderFilters                   variant_filter_type;
	typedef algorithm::Interval<U32, S64>          interval_type;
	typedef io::BasicReader                        basic_reader_type;
	typedef encryption::EncryptionDecorator        encryption_manager_type;

	// Function pointers
	typedef void (self_type::*print_format_function)(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const;
	typedef void (self_type::*print_info_function)(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const;
	typedef void (self_type::*print_filter_function)(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const;
	typedef bool (self_type::*filter_intervals_function)(const meta_entry_type& meta_entry) const;

	typedef buffer_type& (*print_meta_function)(buffer_type& buffer, const char& delimiter, const meta_entry_type& meta_entry, const header_type& header, const block_settings_type& controller);

public:
	VariantReader();
	VariantReader(const std::string& filename);
	VariantReader(const self_type& other);
	virtual ~VariantReader();

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

	/**<
	 * Retrieve the current filter settings for the variant reader. This
	 * object controls the pointers to filter applied to each variant.
	 * @return A reference instance of the filter object
	 */
	inline variant_filter_type& getFilterSettings(void){ return(this->variant_filters); }

	// Basic accessors
	inline header_type& getGlobalHeader(void){ return(this->global_header); }
	inline const header_type& getGlobalHeader(void) const{ return(this->global_header); }
	inline footer_type& getGlobalFooter(void){ return(this->global_footer); }
	inline const footer_type& getGlobalFooter(void) const{ return(this->global_footer); }
	inline index_type& getIndex(void){ return(this->index); }
	inline const index_type& getIndex(void) const{ return(this->index); }
	inline size_t getFilesize(void) const{ return(this->basic_reader.filesize_); }
	inline variant_container_type& getCurrentBlock(void){ return(this->variant_container); }
	inline const variant_container_type& getCurrentBlock(void) const{ return(this->variant_container); }

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
		this->basic_reader.filename_ = filename;
		this->settings.input = filename;
		if(settings.keychain_file.size()){
			if(this->loadKeychainFile() == false)
				return false;
		}
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
	 * Not implemented
	 * @param position
	 * @return
	 */
	bool seektoBlock(const U32 position);

	/**<
	 * Not implemented
	 * @param chromosome_name
	 * @return
	 */
	bool seekToBlockChromosome(const std::string& chromosome_name);

	/**<
	 * Not implemented
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
	bool NextBlock(void);

	/**<
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool GetBlock(const index_entry_type& index_entry);


	/**<
	 * Seeks to a specific YON block without loading anything.
	 * This allows the user to seek to a specific block and
	 * change the block_settings (i.e. what fields to load) and
	 * then invoke nextBlock() for example.
	 * @param blockID
	 * @return
	 */
	bool seek_to_block(const U32& blockID);

	U64 OutputVcf(void);
	U64 OutputVcfLinear(void);
	U64 OutputVcfSearch(void);
	void OuputVcfWrapper(io::BasicBuffer& output_buffer, const yon1_t& entry) const;
	void OutputInfoVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;
	void OutputFormatVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;
	void OutputFilterVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const;

	/**<
	 * Wrapper function to call internal functions `outputCustom` or `outputBlockVCF`.
	 * Decides internally what function to invoke.
	 * @return
	 */
	U64 outputVCF(void);

	// Filter interval intersection and dummy version
	inline bool filterIntervalsDummy(const meta_entry_type& meta_entry) const{ return true; }
	inline bool filterIntervals(const meta_entry_type& meta_entry) const{ return(this->interval_container.find_overlaps(meta_entry).size()); }

	// Calculations
	TACHYON_VARIANT_CLASSIFICATION_TYPE ClassifyVariant(const meta_entry_type& meta, const U32& allele) const;

	/**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	inline bool addIntervals(std::vector<std::string>& interval_strings){
		return(this->interval_container.parseIntervals(interval_strings, this->global_header, this->index));
	}


	/**<
	 *
	 * @param path
	 * @return
	 */
	bool loadKeychainFile(void){
		std::ifstream keychain_reader(settings.keychain_file, std::ios::binary | std::ios::in);
		if(!keychain_reader.good()){
			std::cerr << tachyon::utility::timestamp("ERROR") <<  "Failed to open keychain: " << settings.keychain_file << "..." << std::endl;
			return false;
		}

		keychain_reader >> this->keychain;
		if(!keychain_reader.good()){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse keychain..." << std::endl;
			return false;
		}
		return true;
	}

	/**<
	 *
	 * @param stream
	 */
	void printHeaderVCF(std::ostream& stream = std::cout){
		this->global_header.literals_ += "##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		this->global_header.literals_ += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
		                             +   SSLeay_version(SSLEAY_VERSION) + ","
		                             +  "ZSTD-" + ZSTD_versionString()
		                             +  "; timestamp=" + tachyon::utility::datetime() + "\n";

		this->global_header.literals_ += "##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + "\n";
		this->global_header.literals_ += this->getSettings().get_settings_string();
		this->global_header.literals_ += '\n';

		this->global_header.PrintVcfHeader(stream);
	}


	//<----------------- EXAMPLE FUNCTIONS -------------------------->


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
		if((this->block_settings.load_static & YON_BLK_BV_ALLELES) == false || (this->block_settings.load_static & YON_BLK_BV_GT_INT8) == false || (this->block_settings.load_static & YON_BLK_BV_CONTROLLER) == false || (this->block_settings.load_static & YON_BLK_BV_ID_INFO) == false){
			std::cerr << utility::timestamp("ERROR") << "Cannot run function without loading: SET-MEMBERSHIP, GT, REF or ALT, CONTIG or POSITION..." << std::endl;
			return;
		}

		if(buffer.back() != ';' && buffer.back() != '\t') buffer += ';';

		//objects_type objects;
		//this->loadObjects(objects);
		//U32 n_variants_parsed = 0;

		//for(U32 i = 0; i < objects.genotypes->size(); ++i){
			if(objects.meta_container->at(position).IsDiploid() == false){
				std::cerr << "is not diploid" << std::endl;
				return;
			}

			// If set membership is -1 then calculate all fields
			// Set target FLAG set to all ones; update with actual values if they exist
			U16 target_flag_set = 65535;
			//if(objects.meta_container->at(position).GetInfoPatternId() != -1)
			//	target_flag_set = objects.additional_info_execute_flag_set[objects.meta_container->at(position).GetInfoPatternId()];

			// Get genotype summary data
			objects.genotype_container->at(position).getSummary(*objects.genotype_summary);
			std::vector<double> hwe_p = objects.genotype_summary->calculateHardyWeinberg(objects.meta_container->at(position));
			std::vector<double> af    = objects.genotype_summary->calculateAlleleFrequency(objects.meta_container->at(position));


			//utility::to_vcf_string(stream, this->block_settings.custom_delimiter_char, meta, this->global_header);

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
				buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->ClassifyVariant(objects.meta_container->at(position), 1)];

				for(U32 p = 2; p < objects.meta_container->at(position).n_alleles; ++p){
					buffer += ',';
					buffer += TACHYON_VARIANT_CLASSIFICATION_STRING[this->ClassifyVariant(objects.meta_container->at(position), p)];
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

private:
	basic_reader_type       basic_reader;
	variant_container_type  variant_container;
	block_settings_type     block_settings;
	settings_type           settings;
	variant_filter_type     variant_filters;
	header_type             global_header;
	footer_type             global_footer;
	index_type              index;
	checksum_type           checksums;
	codec_manager_type      codec_manager;
	keychain_type           keychain;
	interval_container_type interval_container;
};

}

#endif /* CORE_TACHYON_READER_H_ */
