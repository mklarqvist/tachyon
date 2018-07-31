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
	bool nextBlock(void);

	/**<
	 * Get the target YON block
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool getBlock(const index_entry_type& index_entry);

	/**<
	 * Get the current YON block in-order as a copy
	 * @return Returns a YON block container. The container has a size of 0 upon fail/empty
	 */
	variant_container_type getBlock(void);


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
	 * Wrapper function to call internal functions `outputCustom` or `outputBlockVCF`.
	 * Decides internally what function to invoke.
	 * @return
	 */
	U64 outputVCF(void);

	/**<
	 *
	 * @return
	 */
	U64 outputCustom(void);

	/**<
	 *
	 * @return
	 */
	U32 outputBlockVCF(void);

	/**<
	 *
	 * @return
	 */
	U32 outputBlockCustom(void);

	// Dummy functions as interfaces for function pointers
	inline void printFILTERDummy(buffer_type& outputBuffer, const U32& position, const objects_type& objects) const{}
	inline void printFORMATDummy(buffer_type& buffer, const char& delimiter, const U32& position, const objects_type& objects, std::vector<core::GTObject>& genotypes_unpermuted) const{}
	inline void printINFODummy(buffer_type& outputBuffer, const char& delimiter, const U32& position, const objects_type& objects) const{}

	// Filter interval intersection and dummy version
	inline bool filterIntervalsDummy(const meta_entry_type& meta_entry) const{ return true; }
	inline bool filterIntervals(const meta_entry_type& meta_entry) const{ return(this->interval_container.find_overlaps(meta_entry).size()); }

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
	 * Not implemented
	 * Return bit-mask primitive of variant classifications detected
	 * @param meta Input meta entry for a site
	 * @return     Returns a primitive interpreted as a boolean presence/absence bit-mask
	 */
	BYTE classifyVariant(const meta_entry_type& meta) const;

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
		this->global_header.literals_ += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		this->global_header.literals_ += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
		                             +   SSLeay_version(SSLEAY_VERSION) + ","
		                             +  "ZSTD-" + ZSTD_versionString()
		                             +  "; timestamp=" + tachyon::utility::datetime();

		this->global_header.literals_ += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + "\n";
		this->global_header.literals_ += this->getSettings().get_settings_string();

		stream << this->global_header.literals_ << std::endl;
		//this->global_header.writeHeaderVCF(stream, true);
	}

	//<----------------- EXAMPLE FUNCTIONS -------------------------->


	U64 timings_meta(){
		containers::MetaContainer meta(this->variant_container.getBlock());
		buffer_type temp(meta.size() * 1000);

		for(U32 p = 0; p < meta.size(); ++p){
			utility::to_vcf_string(temp, this->block_settings.custom_delimiter_char, meta[p], this->global_header);
			//utility::to_vcf_string(std::cout, meta[p], this->global_header);
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
		containers::MetaContainer meta(this->variant_container.getBlock());
		containers::GenotypeContainer gt(this->variant_container.getBlock(), meta);

		for(U32 i = 0; i < gt.size(); ++i){
			// All of these functions are in relative terms very expensive!
			// Avoid using them unless you absolutely have to!
			// Vector of literal genotype representations (lower level)
			//std::vector<core::GTObject> objects     = gt[i].getLiteralObjects();
			// Vector of genotype objects (high level permuted)
			//std::vector<core::GTObject> objects_all = gt[i].getObjects(this->global_header.getSampleNumber());
			// Vector of genotype objects (high level unpermuted - original)
			std::vector<core::GTObject> objects_true = gt[i].getObjects(this->global_header.GetNumberSamples(), this->variant_container.getBlock().ppa_manager);

			std::cout << (int)objects_true[i].alleles[0].allele << (objects_true[i].alleles[1].phase ? '/' : '|') << (int)objects_true[i].alleles[1].allele;
			for(U32 j = 1; j < objects_true.size(); ++j){
				std::cout << '\t' << (int)objects_true[j].alleles[0].allele << (objects_true[j].alleles[1].phase ? '/' : '|') << (int)objects_true[j].alleles[1].allele;
			}
			std::cout << std::endl;
		}
		return(gt.size());
	}

	U64 calculateIBS(math::SquareMatrix<double>& square, math::SquareMatrix<double>& square_temporary){
		algorithm::Timer timer;
		timer.Start();

		containers::MetaContainer meta(this->variant_container.getBlock());
		containers::GenotypeContainer gt(this->variant_container.getBlock(), meta);
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].comparePairwise(square_temporary);

		//square /= (U64)2*this->global_header.getSampleNumber()*gt.size();
		square.addUpperTriagonal(square_temporary, this->variant_container.getBlock().ppa_manager);
		square_temporary.clear();

		// 2 * (Upper triagonal + diagonal) * number of variants
		const U64 updates = 2*((this->global_header.GetNumberSamples()*this->global_header.GetNumberSamples() - this->global_header.GetNumberSamples())/2 + this->global_header.GetNumberSamples()) * gt.size();
		std::cerr << utility::timestamp("DEBUG") << "Updates: " << utility::ToPrettyString(updates) << '\t' << timer.ElapsedString() << '\t' << utility::ToPrettyString((U64)((double)updates/timer.Elapsed().count())) << "/s" << std::endl;
		return((U64)2*this->global_header.GetNumberSamples()*gt.size());
	}

	U64 getTiTVRatios(std::vector<core::TsTvObject>& global){
		containers::MetaContainer meta(this->variant_container.getBlock());
		containers::GenotypeContainer gt(this->variant_container.getBlock(), meta);

		std::vector<core::TsTvObject> objects(this->global_header.GetNumberSamples());
		for(U32 i = 0; i < gt.size(); ++i)
			gt[i].getTsTv(objects);

		for(U32 i = 0; i < objects.size(); ++i)
			global[this->variant_container.getBlock().ppa_manager[i]] += objects[i];

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
			if(objects.meta_container->at(position).IsDiploid() == false){
				std::cerr << "is not diploid" << std::endl;
				return;
			}

			// If set membership is -1 then calculate all fields
			// Set target FLAG set to all ones; update with actual values if they exist
			U16 target_flag_set = 65535;
			if(objects.meta_container->at(position).GetInfoPatternId() != -1)
				target_flag_set = objects.additional_info_execute_flag_set[objects.meta_container->at(position).GetInfoPatternId()];

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

	U64 countVariants(void){
		containers::MetaContainer meta(this->variant_container.getBlock());
		return(meta.size());
	}

	U64 iterateMeta(std::ostream& stream = std::cout){
		containers::MetaContainer meta(this->variant_container.getBlock());
		containers::GenotypeContainer gt(this->variant_container.getBlock(), meta);
		containers::GenotypeSummary gt_summary;
		for(U32 i = 0; i < gt.size(); ++i){
			// If there's > 5 alleles continue
			if(gt[i].getMeta().GetNumberAlleles() >= 5) continue;
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

		const io::VcfInfo* info = this->global_header.GetInfo("AF");
		if(info != nullptr){
			containers::InfoContainer<double> it_i(this->variant_container.getBlock().info_containers[1]);
			//math::MathSummaryStatistics stats = it_i.getSummaryStatistics();
			//std::cerr << stats.n_total << '\t' << stats.mean << '\t' << stats.standard_deviation << '\t' << stats.min << "-" << stats.max << std::endl;
			for(U32 i = 0; i < it_i.size(); ++i){
				//if(it_i[i].size() < 3) continue;
				//it[i].toVCFString(stream, this->global_header, this->variant_container.getBlock().index_entry.contigID, this->variant_container.getBlock().index_entry.minPosition);

				//stream << (int)it_i[i][0];
				for(U32 j = 0; j < it_i[i].size(); ++j)
					stream << it_i[i][j] << ' ';
			}
			stream << '\n';
		}
		return(0);
	}

	/**<
	 * Construct Occ matrix from a target input file
	 * @param file
	 * @return
	 */
	bool loadGroups(const std::string& file){
		/*
		if(this->Occ.size() != 0){
			std::cerr << "groups already loaded" << std::endl;
			return false;
		}

		if(file.size() == 0){
			std::cerr << "No file set" << std::endl;
			return false;
		}

		// Open stream
		std::ifstream stream(file);
		if(!stream.good()){
			std::cerr << "bad file" << std::endl;
			return false;
		}

		this->group_htable = new hash_table(10000);


		std::string line;
		std::vector< std::vector< S32 > > groupings(this->samples);

		S32* sampleID = nullptr;
		S32* groupID_lookup = nullptr;
		S32 groupID = 0;
		U32 n_lines = 0;

		while(getline(stream, line)){
			// Empty lines
			if(line.size() == 0)
				break;

			// Assert correct format
			// Count tabs until out-of-range
			size_t prev_pos = 0;
			size_t pos = line.find('\t', prev_pos + 1);

			if(pos == std::string::npos){
				std::cerr << "illegal: has no data for sample (" << n_lines << ")" << std::endl;
				return(false);
			}

			// Make sure the sample name exists in this file
			const std::string sampleName = std::string(&line[prev_pos], pos - prev_pos);
			if(!this->totempole.sampleHashTable->GetItem(&sampleName[0], &sampleName, sampleID, sampleName.length())){
				std::cerr << "sample does not exist" << std::endl;
				return false;
			}

			// Cycle over lines
			prev_pos = pos + 1;
			while(prev_pos != std::string::npos){
				pos = line.find('\t', prev_pos + 1);
				if(pos == std::string::npos){
					pos = line.size();
				}

				const std::string group = std::string(&line[prev_pos], pos - prev_pos);
				if(!this->group_htable->GetItem(&group[0], &group, groupID_lookup, group.length())){
					this->group_htable->SetItem(&group[0], &group, groupID, group.length());
					this->groups.push_back(GroupPair(group));

					groupings[*sampleID].push_back(groupID);
					++groupID;
				} else {
					groupings[*sampleID].push_back(*groupID_lookup);
					++this->groups[*groupID_lookup];
				}

				if(pos == line.size())
					break;

				prev_pos = pos + 1;
			}

			++n_lines;
		}

		if(groupID == 0){
			std::cerr << "no data loaded" << std::endl;
			return false;
		}

		this->Occ = occ_matrix(this->samples + 1, std::vector< U64 >(groupID, 0));
		occ_vector* prev = &Occ[0];

		for(U32 i = 0; i < this->samples; ++i){
			// Propagate previous vector counts
			for(U32 j = 0; j < groupID; ++j)
				this->Occ[i + 1][j] = prev->at(j);

			// Update
			for(U32 j = 0; j < groupings[i].size(); ++j)
				this->Occ[i + 1][groupings[i][j]] = prev->at(groupings[i][j]) + 1;

			prev = &this->Occ[i + 1];
		}

		for(U32 i = 0; i < this->groups.size(); ++i)
			std::cerr << i << '\t' << this->groups[i].name << '\t' << this->groups[i].n_entries << std::endl;

		// Temp
		// Dump
		 */
		/*
		for(U32 i = 0; i < this->Occ.size(); ++i){
			for(U32 j = 0; j < this->Occ[0].size(); ++j){
				std::cout << this->Occ[i][j] << '\t';
			}
			std::cout << std::endl;
		}
		*/

		return(true);
	}

protected:
	basic_reader_type       basic_reader;
	//std::ifstream           stream;
	//size_t                  filesize;
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
