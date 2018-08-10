#include "variant_reader.h"

namespace tachyon{

VariantReader::VariantReader()
{}

VariantReader::VariantReader(const std::string& filename) :
	basic_reader(filename)
{}

VariantReader::~VariantReader(){}

VariantReader::VariantReader(const self_type& other) :
	basic_reader(other.basic_reader),
	block_settings(other.block_settings),
	settings(other.settings),
	global_header(other.global_header),
	global_footer(other.global_footer),
	index(other.index),
	checksums(other.checksums),
	keychain(other.keychain)
{
	this->basic_reader.open();
}

bool VariantReader::open(void){
	if(this->settings.input.size() == 0){
		std::cerr << utility::timestamp("ERROR") << "No input file specified!" << std::endl;
		return false;
	}

	if(this->basic_reader.open() == false){
		std::cerr << "Failed to open" << std::endl;
		return false;
	}

	if(this->basic_reader.filesize_ <= YON_FOOTER_LENGTH){
		std::cerr << utility::timestamp("ERROR") << "File is corrupted!" << std::endl;
		return false;
	}

	// Seek to start of footer
	this->basic_reader.stream_.seekg((U64)this->basic_reader.filesize_ - YON_FOOTER_LENGTH);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to seek in file!" << std::endl;
		return false;
	}
	this->basic_reader.stream_ >> this->global_footer;

	// Validate footer
	if(this->global_footer.validate() == false){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate footer!" << std::endl;
		return false;
	}

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to read file!" << std::endl;
		return false;
	}

	// Seek to start of file
	this->basic_reader.stream_.seekg(0);
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to rewind file!" << std::endl;
		return false;
	}

	// Load header
	//this->stream >> this->global_header;
	char magic_string[tachyon::constants::FILE_HEADER_LENGTH];
	this->basic_reader.stream_.read(&magic_string[0], tachyon::constants::FILE_HEADER_LENGTH);
	if(strncmp(&magic_string[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
		std::cerr << utility::timestamp("ERROR") << "Failed to validate Tachyon magic string!" << std::endl;
		return false;
	}

	uint32_t l_data   = 0;
	uint32_t l_c_data = 0;
	utility::DeserializePrimitive(l_data, this->basic_reader.stream_);
	utility::DeserializePrimitive(l_c_data, this->basic_reader.stream_);

	io::BasicBuffer header_uncompressed(l_data + 1024);
	io::BasicBuffer header_compressed(l_c_data + 1024); header_compressed.n_chars   = l_c_data;

	this->basic_reader.stream_.read(header_compressed.data(), l_c_data);

	if(!this->codec_manager.zstd_codec.Decompress(header_compressed, header_uncompressed)){
		std::cerr << utility::timestamp("ERROR") << "Failed to decompress header!" << std::endl;
		return false;
	}
	assert(header_uncompressed.size() == l_data);
	header_uncompressed >> this->global_header; // parse header from buffer

	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR") << "Failed to get header!" << std::endl;
		return false;
	}

	this->variant_container << this->global_header;

	// Keep track of start position
	const U64 return_pos = this->basic_reader.stream_.tellg();
	this->basic_reader.stream_.seekg(this->global_footer.offset_end_of_data);
	this->basic_reader.stream_ >> this->index;
	this->basic_reader.stream_ >> this->checksums;
	this->basic_reader.stream_.seekg(return_pos);

	return(this->basic_reader.stream_.good());
}

bool VariantReader::NextBlock(){
	// If the stream is faulty then return
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// If the current position is the EOF then
	// exit the function
	if((U64)this->basic_reader.stream_.tellg() == this->global_footer.offset_end_of_data)
		return false;

	// Reset and re-use
	this->variant_container.reset();

	if(!this->variant_container.GetBlock().ReadHeaderFooter(this->basic_reader.stream_))
		return false;


	if(!this->codec_manager.zstd_codec.Decompress(this->variant_container.GetBlock().footer_support)){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression of footer!" << std::endl;
	}
	this->variant_container.GetBlock().footer_support.buffer_data_uncompressed >> this->variant_container.GetBlock().footer;

	// Attempts to read a YON block with the settings provided
	if(!this->variant_container.ReadBlock(this->basic_reader.stream_, this->block_settings))
		return false;

	// encryption manager ascertainment
	if(this->variant_container.AnyEncrypted()){
		if(this->keychain.size() == 0){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Data is encrypted but no keychain was provided!" << std::endl;
			return false;
		}

		encryption_manager_type encryption_manager;
		if(!encryption_manager.decryptAES256(this->variant_container.GetBlock(), this->keychain)){
			std::cerr << utility::timestamp("ERROR", "DECRYPTION") << "Failed decryption!" << std::endl;
			return false;
		}
	}

	// Internally decompress available data
	if(!this->codec_manager.Decompress(this->variant_container.GetBlock())){
		std::cerr << utility::timestamp("ERROR", "COMPRESSION") << "Failed decompression!" << std::endl;
		return false;
	}

	// All passed
	return true;
}

bool VariantReader::GetBlock(const index_entry_type& index_entry){
	// If the stream is not good then return.
	if(!this->basic_reader.stream_.good()){
		std::cerr << utility::timestamp("ERROR", "IO") << "Corrupted! Input stream died prematurely!" << std::endl;
		return false;
	}

	// Seek to target block id with the help of the linear index.
	this->basic_reader.stream_.seekg(index_entry.byte_offset);
	if(this->basic_reader.stream_.good() == false){
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to seek to given offset using target index entry!" << std::endl;
		return(false);
	}

	// Load the next block if possible.
	return(this->NextBlock());
}

TACHYON_VARIANT_CLASSIFICATION_TYPE VariantReader::ClassifyVariant(const meta_entry_type& meta, const U32& allele) const{
	const S32 ref_size = meta.alleles[0].size();
	const S32 diff = ref_size - meta.alleles[allele].size();
	//std::cerr << diff << ",";
	if(meta.alleles[0].allele[0] == '<' || meta.alleles[allele].allele[0] == '<') return(YON_VARIANT_CLASS_SV);
	else if(diff == 0){
		if(ref_size == 1 && meta.alleles[0].allele[0] != meta.alleles[allele].allele[0]){
			if(meta.alleles[allele].allele[0] == 'A' || meta.alleles[allele].allele[0] == 'T' || meta.alleles[allele].allele[0] == 'G' || meta.alleles[allele].allele[0] == 'C')
				return(YON_VARIANT_CLASS_SNP);
			else return(YON_VARIANT_CLASS_UNKNOWN);
		}
		else if(ref_size != 1){
			U32 characters_identical = 0;
			const U32 length_shortest = ref_size < meta.alleles[allele].size() ? ref_size : meta.alleles[allele].size();

			for(U32 c = 0; c < length_shortest; ++c){
				characters_identical += (meta.alleles[0].allele[c] == meta.alleles[allele].allele[c]);
			}

			if(characters_identical == 0) return(YON_VARIANT_CLASS_MNP);
			else return(YON_VARIANT_CLASS_CLUMPED);
		}
	} else {
		const U32 length_shortest = ref_size < meta.alleles[allele].size() ? ref_size : meta.alleles[allele].size();
		U32 characters_non_standard = 0;
		for(U32 c = 0; c < length_shortest; ++c){
			characters_non_standard += (meta.alleles[allele].allele[c] != 'A' && meta.alleles[allele].allele[c] != 'T' && meta.alleles[allele].allele[c] != 'C' && meta.alleles[allele].allele[c] !='G');
		}
		if(characters_non_standard) return(YON_VARIANT_CLASS_UNKNOWN);
		else return(YON_VARIANT_CLASS_INDEL);
	}
	return(YON_VARIANT_CLASS_UNKNOWN);
}

void VariantReader::OuputVcfWrapper(io::BasicBuffer& output_buffer, yon1_t& entry) const{
	utility::to_vcf_string(output_buffer, '\t', *entry.meta, this->global_header);
	output_buffer += '\t';

	// Print Filter, Info, and Format if available.
	this->OutputFilterVcf(output_buffer, entry);
	this->OutputInfoVcf(output_buffer, entry);
	this->OutputFormatVcf(output_buffer, entry);
	output_buffer += '\n';

	if(output_buffer.size() > 65536){
		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
	}
}

void VariantReader::OutputInfoVcf(io::BasicBuffer& output_buffer, yon1_t& entry) const{
	// Print Info.
	if(entry.n_info){
		const uint32_t n_info_avail = entry.info_ids->size();
		if(n_info_avail){
			if(entry.info_hdr[0]->yon_type == YON_VCF_HEADER_FLAG){
				output_buffer += entry.info_hdr[0]->id;
			} else {
				output_buffer += entry.info_hdr[0]->id;
				output_buffer += '=';
				entry.info[0]->to_vcf_string(output_buffer);
			}

			for(U32 j = 1; j < n_info_avail; ++j){
				output_buffer += ';';
				if(entry.info_hdr[j]->yon_type == YON_VCF_HEADER_FLAG){
					output_buffer += entry.info_hdr[j]->id;
				} else {
					output_buffer += entry.info_hdr[j]->id;
					output_buffer += '=';
					entry.info[j]->to_vcf_string(output_buffer);
				}
			}

			if(this->GetBlockSettings().annotate_extra){
				entry.EvaluateSummary(true);
				entry.gt_sum->d->PrintVcf(output_buffer);
			}
		}
	} else {
		if(this->GetBlockSettings().annotate_extra){
			entry.EvaluateSummary(true);
			entry.gt_sum->d->PrintVcf(output_buffer);
		} else
			output_buffer += '.';
	}
}

void VariantReader::OutputFormatVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const{
	if(entry.n_format){
		output_buffer += '\t';

		const uint32_t n_format_avail = entry.format_ids->size();
		if(n_format_avail){
			output_buffer += entry.format_hdr[0]->id;
			for(U32 j = 1; j < n_format_avail; ++j){
				output_buffer += ':';
				output_buffer += entry.format_hdr[j]->id;
			}
			output_buffer += '\t';

			// Vcf FORMAT values are interleaved such that values are
			// presented in a columnar representation with data for
			// each sample is concatenated together such as the pattern
			// GT:AF:AS is display for three samples as:
			//
			// Sample 1    Sample 2    Sample 3
			// 0|0:0.5|ABC 0|0:0.5|ABC 0|0:0.5|ABC
			//
			// This memory layout requires the interleaving of the
			// internally separated data streams resulting in slower
			// Vcf printing speeds compared to the naive Bcf file format.
			//
			// First calculate the FORMAT:GT field for this variant site.
			// Case when the only available FORMAT field is the GT field.
			if(n_format_avail == 1 && entry.is_loaded_gt && entry.meta->controller.gt_available){
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				// Iterate over samples and print FORMAT:GT value in Vcf format.
				entry.gt->d_exp[0]->PrintVcf(output_buffer, entry.gt->m);
				for(U32 s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.gt->d_exp[s]->PrintVcf(output_buffer, entry.gt->m);
				}

				entry.gt->d_exp = nullptr;
			}
			// Case when there are > 1 Vcf Format fields and the GT field
			// is available.
			else if(n_format_avail > 1 && entry.is_loaded_gt && entry.meta->controller.gt_available){
				entry.gt->ExpandExternal(this->variant_container.GetAllocatedGenotypeMemory());
				entry.gt->d_exp = this->variant_container.GetAllocatedGenotypeMemory();

				entry.gt->d_exp[0]->PrintVcf(output_buffer, entry.gt->m);
				for(U32 g = 1; g < n_format_avail; ++g){
					output_buffer += ':';
					entry.fmt[g]->to_vcf_string(output_buffer, 0);
				}
				for(U32 s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.gt->d_exp[s]->PrintVcf(output_buffer, entry.gt->m);
					for(U32 g = 1; g < n_format_avail; ++g){
						output_buffer += ':';
						entry.fmt[g]->to_vcf_string(output_buffer, s);
					}
				}

				entry.gt->d_exp = nullptr;
			}
			// All other cases.
			else {
				entry.fmt[0]->to_vcf_string(output_buffer, 0);
				for(U32 g = 1; g < n_format_avail; ++g){
					output_buffer += ':';
					entry.fmt[g]->to_vcf_string(output_buffer, 0);
				}

				for(U32 s = 1; s < this->global_header.GetNumberSamples(); ++s){
					output_buffer += '\t';
					entry.fmt[0]->to_vcf_string(output_buffer, s);
					for(U32 g = 1; g < n_format_avail; ++g){
						output_buffer += ':';
						entry.fmt[g]->to_vcf_string(output_buffer, s);
					}
				}
			}
		}
	}
}

void VariantReader::OutputFilterVcf(io::BasicBuffer& output_buffer, const yon1_t& entry) const{
	if(entry.n_filter){
		const uint32_t n_filter_avail = entry.filter_ids->size();
		if(n_filter_avail){
			output_buffer += entry.filter_hdr[0]->id;
			for(U32 j = 1; j < n_filter_avail; ++j){
				output_buffer += ';';
				output_buffer += entry.filter_hdr[j]->id;
			}
		} else {
			output_buffer += '.';
		}
	} else output_buffer += '.';
	output_buffer += '\t';
}

U64 VariantReader::OutputVcfLinear(void){

	this->variant_container.AllocateGenotypeMemory();

	while(this->NextBlock()){
		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);
		io::BasicBuffer output_buffer(100000);

		for(U32 i = 0; i < objects->meta_container->size(); ++i){
			if(this->variant_filters.filter(entries[i], i) == false)
				continue;

			this->OuputVcfWrapper(output_buffer, entries[i]);
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		delete [] entries;
		delete objects;
	}
	return 0;
}

U64 VariantReader::OutputVcfSearch(void){
	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervals;

	for(U32 i = 0; i < this->interval_container.getBlockList().size(); ++i){
		this->GetBlock(this->interval_container.getBlockList()[i]);

		objects_type* objects = this->GetCurrentContainer().LoadObjects(this->block_settings);
		yon1_t* entries = this->GetCurrentContainer().LazyEvaluate(*objects);
		io::BasicBuffer output_buffer(100000);

		for(U32 i = 0; i < objects->meta_container->size(); ++i){
			if((this->*filter_intervals)(objects->meta_container->at(i)) == false)
				continue;

			if(this->variant_filters.filter(entries[i], i) == false)
				continue;

			this->OuputVcfWrapper(output_buffer, entries[i]);
		}

		std::cout.write(output_buffer.data(), output_buffer.size());
		output_buffer.reset();
		delete [] entries;
		delete objects;
	}

	return 0;
}

U64 VariantReader::OutputVcf(void){
	if(this->GetBlockSettings().show_vcf_header)
		this->global_header.PrintVcfHeader(std::cout);

	this->interval_container.build(this->global_header);

	// Filter functionality
	filter_intervals_function filter_intervals = &self_type::FilterIntervalsDummy;
	if(this->interval_container.size()) return(this->OutputVcfSearch());
	else return(this->OutputVcfLinear());
}

/**<
 * Outputs
 * @return
 */
U64 VariantReader::outputVCF(void){
	U64 n_variants = 0;

	if(this->block_settings.annotate_extra){
		// fixme
		// if special
		// "FS_A", "AN", "NM", "NPM", "AC", "AC_FW", "AC_REV", "AF", "HWE_P", "VT", "MULTI_ALLELIC", "F_PIC"
		if(this->global_header.GetInfo("FS_A") == nullptr)          this->global_header.literals_ += "\n##INFO=<ID=FS_A,Number=A,Type=Float>";
		if(this->global_header.GetInfo("AN") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AN,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("NM") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=NM,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("NPM") == nullptr)           this->global_header.literals_ += "\n##INFO=<ID=NPM,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AC,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC_FWD") == nullptr)        this->global_header.literals_ += "\n##INFO=<ID=AC_FWD,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("AC_REV") == nullptr)        this->global_header.literals_ += "\n##INFO=<ID=AC_REV,Number=A,Type=Integer>";
		if(this->global_header.GetInfo("HWE_P") == nullptr)         this->global_header.literals_ += "\n##INFO=<ID=HWE_P,Number=A,Type=Float>";
		if(this->global_header.GetInfo("VT") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=VT,Number=A,Type=String>";
		if(this->global_header.GetInfo("AF") == nullptr)            this->global_header.literals_ += "\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"Estimated allele frequency in the range (0,1)\">";
		if(this->global_header.GetInfo("MULTI_ALLELIC") == nullptr) this->global_header.literals_ += "\n##INFO=<ID=MULTI_ALLELIC,Number=0,Type=Flag>";
		if(this->global_header.GetInfo("F_PIC") == nullptr)         this->global_header.literals_ += "\n##INFO=<ID=F_PIC,Number=A,Type=Float,Description=\"Population inbreeding coefficient (F-statistics)\">";
	}

	this->global_header.literals_ += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
	this->global_header.literals_ += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
			  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + utility::datetime();

	this->global_header.literals_ += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + '\n';
	this->global_header.literals_ += this->GetSettings().get_settings_string();

	// Output VCF header
	if(this->block_settings.show_vcf_header){
		//this->global_header.writeHeaderVCF(std::cout, this->block_settings.format_all.load || this->block_settings.format_list.size());
	}

	// If seek is active for targetted intervals
	if(this->interval_container.hasIntervals()){
		if(this->interval_container.build(this->global_header) == false)
			return false;

		if(this->interval_container.getBlockList().size()){
			for(U32 i = 0; i < this->interval_container.getBlockList().size(); ++i){
				if(this->GetBlock(this->interval_container.getBlockList()[i]) == false){
					return(0);
				}
				//n_variants += this->outputBlockVCF();
			}

			return(n_variants);
		} else { // Provided intervals but no matching YON blocks
			return(0);
		}
	}

	// While there are YON blocks
	//while(this->NextBlock()) n_variants += this->outputBlockVCF();
	return(n_variants);
}

}
