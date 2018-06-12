#include "variantblock.h"

#include "../support/helpers.h"
#include "../support/type_definitions.h"
#include "../algorithm/compression/compression_container.h"

namespace tachyon{
namespace containers{

VariantBlock::VariantBlock() :
	info_containers(new container_type[200]),
	format_containers(new container_type[200]),
	end_block_(0),
	start_compressed_data_(0),
	end_compressed_data_(0),
	n_info_loaded(0),
	n_format_loaded(0)
{
	// Base container streams are always of type TYPE_STRUCT
	this->meta_alleles_container.setType(YON_TYPE_STRUCT);
	this->meta_controller_container.setType(tachyon::YON_TYPE_16B);
	this->meta_refalt_container.setType(tachyon::YON_TYPE_8B);

	this->footer_support.resize(65536);
}

VariantBlock::~VariantBlock(){
	delete [] this->info_containers;
	delete [] this->format_containers;
}

void VariantBlock::clear(void){
	for(U32 i = 0; i < this->footer.n_info_streams; ++i)
		this->info_containers[i].reset();

	for(U32 i = 0; i < this->footer.n_format_streams; ++i)
		this->format_containers[i].reset();

	this->header.reset();
	this->footer.reset();
	this->footer_support.reset();

	this->meta_contig_container.reset();
	this->meta_positions_container.reset();
	this->meta_refalt_container.reset();
	this->meta_controller_container.reset();
	this->meta_quality_container.reset();
	this->meta_names_container.reset();
	this->meta_alleles_container.reset();
	this->meta_info_map_ids.reset();
	this->meta_format_map_ids.reset();
	this->meta_filter_map_ids.reset();
	this->gt_support_data_container.reset();
	this->gt_rle8_container.reset();
	this->gt_rle16_container.reset();
	this->gt_rle32_container.reset();
	this->gt_rle64_container.reset();
	this->gt_simple8_container.reset();
	this->gt_simple16_container.reset();
	this->gt_simple32_container.reset();
	this->gt_simple64_container.reset();

	// Base container data types are always TYPE_STRUCT
	// Map ID fields are always S32 fields
	this->meta_alleles_container.setType(YON_TYPE_STRUCT);
	this->meta_controller_container.setType(tachyon::YON_TYPE_16B);
	this->meta_refalt_container.setType(tachyon::YON_TYPE_8B);

	this->info_fields.clear();
	this->format_fields.clear();
	this->filter_fields.clear();
	this->info_patterns.clear();
	this->format_patterns.clear();
	this->filter_patterns.clear();

	this->end_block_        = 0;
	this->n_info_loaded       = 0;
	this->n_format_loaded     = 0;
	this->start_compressed_data_       = 0;
	this->end_compressed_data_ = 0;

	this->ppa_manager.reset();
}

void VariantBlock::resize(const U32 s){
	if(s == 0) return;

	this->meta_contig_container.resize(s);
	this->meta_positions_container.resize(s);
	this->meta_refalt_container.resize(s);
	this->meta_controller_container.resize(s);
	this->meta_quality_container.resize(s);
	this->meta_names_container.resize(s);
	this->meta_alleles_container.resize(s);
	this->meta_info_map_ids.resize(s);
	this->meta_format_map_ids.resize(s);
	this->meta_filter_map_ids.resize(s);
	this->gt_support_data_container.resize(s);
	this->gt_rle8_container.resize(s);
	this->gt_rle16_container.resize(s);
	this->gt_rle32_container.resize(s);
	this->gt_rle64_container.resize(s);
	this->gt_simple8_container.resize(s);
	this->gt_simple16_container.resize(s);
	this->gt_simple32_container.resize(s);
	this->gt_simple64_container.resize(s);

	for(U32 i = 0; i < 200; ++i){
		this->info_containers[i].resize(s);
		this->format_containers[i].resize(s);
	}
}

void VariantBlock::updateContainers(void){
	this->meta_contig_container.updateContainer();
	this->meta_positions_container.updateContainer();
	this->meta_refalt_container.updateContainer();
	this->meta_controller_container.updateContainer(false);
	this->meta_quality_container.updateContainer();
	this->meta_names_container.updateContainer();
	this->meta_alleles_container.updateContainer();
	this->meta_filter_map_ids.updateContainer();
	this->meta_format_map_ids.updateContainer();
	this->meta_info_map_ids.updateContainer();
	this->gt_support_data_container.updateContainer();
	this->gt_rle8_container.updateContainer(false);
	this->gt_rle16_container.updateContainer(false);
	this->gt_rle32_container.updateContainer(false);
	this->gt_rle64_container.updateContainer(false);
	this->gt_simple8_container.updateContainer(false);
	this->gt_simple16_container.updateContainer(false);
	this->gt_simple32_container.updateContainer(false);
	this->gt_simple64_container.updateContainer(false);

	for(U32 i = 0; i < this->footer.n_info_streams; ++i){
		assert(this->info_containers[i].header.data_header.stride != 0);
		this->info_containers[i].updateContainer();
	}

	for(U32 i = 0; i < this->footer.n_format_streams; ++i){
		assert(this->format_containers[i].header.data_header.stride != 0);
		this->format_containers[i].updateContainer();
	}
}

bool VariantBlock::readHeaderFooter(std::ifstream& stream){
	if(!stream.good()){
		std::cerr << utility::timestamp("ERROR") << "File stream is corrupted..." << std::endl;
		return false;
	}

	stream >> this->header; // load header
	this->start_compressed_data_ = (U64)stream.tellg(); // start of compressed data
	stream.seekg(this->start_compressed_data_ + this->header.l_offset_footer); // seek to start of footer
	this->end_compressed_data_ = stream.tellg(); // end of compressed data

	U32 footer_uLength = 0;
	U32 footer_cLength = 0;
	U32 footer_crc = 0;
	stream.read(reinterpret_cast<char*>(&footer_uLength), sizeof(U32));
	stream.read(reinterpret_cast<char*>(&footer_cLength), sizeof(U32));
	stream.read(reinterpret_cast<char*>(&footer_crc),     sizeof(U32));
	this->footer_support.resize(footer_cLength);
	stream.read(this->footer_support.buffer_data.data(), footer_cLength);
	this->footer_support.buffer_data.n_chars = footer_cLength;
	this->footer_support.buffer_data_uncompressed.resize(footer_uLength);
	this->footer_support.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	this->footer_support.header.data_header.cLength = footer_cLength;
	this->footer_support.header.data_header.uLength = footer_uLength;
	this->footer_support.header.data_header.crc     = footer_crc;

	// Assert end-of-block marker
	U64 eof_marker;
	stream.read(reinterpret_cast<char*>(&eof_marker), sizeof(U64));
	assert(eof_marker == constants::TACHYON_BLOCK_EOF);
	this->end_block_ = stream.tellg(); // end-of-block offset
	stream.seekg(this->start_compressed_data_);
	return(stream.good());
}

bool VariantBlock::read(std::ifstream& stream, settings_type& settings){
	if(settings.load_ppa){
		if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
			this->ppa_manager.header = this->footer.offset_ppa;
			stream.seekg(this->start_compressed_data_ + this->footer.offset_ppa.data_header.offset);
			stream >> this->ppa_manager;
		}
	}

	if(settings.load_contig){
		this->__loadContainerSeek(stream, this->footer.offset_meta_contig, this->meta_contig_container);
	}

	if(settings.load_positons){
		this->__loadContainerSeek(stream, this->footer.offset_meta_position, this->meta_positions_container);
	}

	if(settings.load_controller){
		this->__loadContainerSeek(stream, this->footer.offset_meta_controllers, this->meta_controller_container);
	}

	if(settings.load_quality){
		this->__loadContainerSeek(stream, this->footer.offset_meta_quality, this->meta_quality_container);
	}

	if(settings.load_names){
		this->__loadContainerSeek(stream, this->footer.offset_meta_names, this->meta_names_container);
	}

	if(settings.load_alleles){
		this->__loadContainerSeek(stream, this->footer.offset_meta_refalt, this->meta_refalt_container);
		this->__loadContainer(stream, this->footer.offset_meta_alleles, this->meta_alleles_container);
	}

	if(settings.load_genotypes_rle){
		this->__loadContainerSeek(stream, this->footer.offset_gt_8b, this->gt_rle8_container);
		this->__loadContainer(stream, this->footer.offset_gt_16b, this->gt_rle16_container);
		this->__loadContainer(stream, this->footer.offset_gt_32b, this->gt_rle32_container);
		this->__loadContainer(stream, this->footer.offset_gt_64b, this->gt_rle64_container);
	}

	if(settings.load_genotypes_simple){
		this->__loadContainerSeek(stream, this->footer.offset_gt_simple8, this->gt_simple8_container);
		this->__loadContainer(stream, this->footer.offset_gt_simple16, this->gt_simple16_container);
		this->__loadContainer(stream, this->footer.offset_gt_simple32, this->gt_simple32_container);
		this->__loadContainer(stream, this->footer.offset_gt_simple64, this->gt_simple64_container);
	}

	if(settings.load_genotypes_support){
		this->__loadContainerSeek(stream, this->footer.offset_gt_helper, this->gt_support_data_container);
	}

	if(settings.load_set_membership){
		this->__loadContainerSeek(stream, this->footer.offset_meta_info_id, this->meta_info_map_ids);
		this->__loadContainer(stream, this->footer.offset_meta_filter_id, this->meta_filter_map_ids);
		this->__loadContainer(stream, this->footer.offset_meta_format_id, this->meta_format_map_ids);
	}

	// Load all info
	if(settings.load_info && this->footer.n_info_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.info_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_info_streams; ++i){
			this->__loadContainer(stream, this->footer.info_offsets[i], this->info_containers[i]);
			++this->n_info_loaded;
			settings.load_info_ID_loaded.push_back(SettingsMap(i,i,&this->footer.info_offsets[i]));
		}
	}
	// If we have supplied a list of identifiers
	else if(settings.load_info_ID_loaded.size()){
		// Ascertain that random access is linearly forward
		for(U32 i = 0; i < settings.load_info_ID_loaded.size(); ++i){
			stream.seekg(this->start_compressed_data_ + settings.load_info_ID_loaded[i].offset->data_header.offset);
			if(!stream.good()){
				std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for INFO field in block!" << std::endl;
				return false;
			}

			//std::cerr << this->footer.info_offsets[settings.load_info_ID_loaded[i].target_stream_local].data_header.uLength << std::endl;

			// Read data
			this->info_containers[settings.load_info_ID_loaded[i].iterator_index].header = this->footer.info_offsets[settings.load_info_ID_loaded[i].target_stream_local];
			stream >> this->info_containers[settings.load_info_ID_loaded[i].iterator_index];
			//std::cerr << "into: " << settings.load_info_ID_loaded[i].iterator_index << std::endl;
			++this->n_info_loaded;
		}
	} // end case load_info_ID

	// Load all FORMAT data
	if(settings.load_format && this->footer.n_format_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.format_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_format_streams; ++i){
			this->__loadContainer(stream, this->footer.format_offsets[i], this->format_containers[i]);
			++this->n_format_loaded;
			settings.load_format_ID_loaded.push_back(SettingsMap(i,i,&this->footer.format_offsets[i]));
		}
		//std::cerr << this->end_compressed_data_ << '/' << (U64)stream.tellg() << std::endl;
		assert(this->end_compressed_data_ == (U64)stream.tellg());
	} // If we have supplied a list of identifiers
	else if(settings.load_format_ID_loaded.size()){
		// Ascertain that random access is linearly forward
		for(U32 i = 0; i < settings.load_format_ID_loaded.size(); ++i){
			stream.seekg(this->start_compressed_data_ + settings.load_format_ID_loaded[i].offset->data_header.offset);
			if(!stream.good()){
				std::cerr << utility::timestamp("ERROR","IO") << "Failed to seek for FORMAT field in block!" << std::endl;
				return false;
			}

			//std::cerr << this->footer.info_offsets[settings.load_info_ID_loaded[i].target_stream_local].data_header.uLength << std::endl;

			// Read data
			this->format_containers[settings.load_format_ID_loaded[i].iterator_index].header = this->footer.format_offsets[settings.load_format_ID_loaded[i].target_stream_local];
			stream >> this->format_containers[settings.load_format_ID_loaded[i].iterator_index];
			//std::cerr << "into: " << settings.load_info_ID_loaded[i].iterator_index << std::endl;
			++this->n_format_loaded;
		}
	} // end case load_info_ID

	stream.seekg(this->end_block_); // seek to end-of-block
	return(true);
}

const U64 VariantBlock::__determineCompressedSize(void) const{
	U64 total = 0;
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		total += this->ppa_manager.getObjectSize();

	total += this->meta_contig_container.getObjectSize();
	total += this->meta_positions_container.getObjectSize();
	total += this->meta_refalt_container.getObjectSize();
	total += this->meta_controller_container.getObjectSize();
	total += this->meta_quality_container.getObjectSize();
	total += this->meta_names_container.getObjectSize();
	total += this->meta_alleles_container.getObjectSize();
	total += this->meta_info_map_ids.getObjectSize();
	total += this->meta_format_map_ids.getObjectSize();
	total += this->meta_filter_map_ids.getObjectSize();
	total += this->gt_support_data_container.getObjectSize();
	total += this->gt_rle8_container.getObjectSize();
	total += this->gt_rle16_container.getObjectSize();
	total += this->gt_rle32_container.getObjectSize();
	total += this->gt_rle64_container.getObjectSize();
	total += this->gt_simple8_container.getObjectSize();
	total += this->gt_simple16_container.getObjectSize();
	total += this->gt_simple32_container.getObjectSize();
	total += this->gt_simple64_container.getObjectSize();

	for(U32 i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].getObjectSize();
	for(U32 i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].getObjectSize();

	return(total);
}

void VariantBlock::updateOutputStatistics(import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format){
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
		stats_basic[1].cost_uncompressed += this->ppa_manager.header.data_header.uLength;
		stats_basic[1].cost_compressed   += this->ppa_manager.header.data_header.cLength;
	}

	stats_basic[2]  += this->meta_contig_container;
	stats_basic[3]  += this->meta_positions_container;
	stats_basic[4]  += this->meta_refalt_container;
	stats_basic[5]  += this->meta_controller_container;
	stats_basic[6]  += this->meta_quality_container;
	stats_basic[7]  += this->meta_names_container;
	stats_basic[8]  += this->meta_alleles_container;
	stats_basic[9]  += this->meta_info_map_ids;
	stats_basic[10] += this->meta_format_map_ids;
	stats_basic[11] += this->meta_filter_map_ids;
	stats_basic[12] += this->gt_support_data_container;
	stats_basic[13] += this->gt_rle8_container;
	stats_basic[14] += this->gt_rle16_container;
	stats_basic[15] += this->gt_rle32_container;
	stats_basic[16] += this->gt_rle64_container;
	stats_basic[17] += this->gt_simple8_container;
	stats_basic[18] += this->gt_simple16_container;
	stats_basic[19] += this->gt_simple32_container;
	stats_basic[20] += this->gt_simple64_container;

	for(U32 i = 0; i < this->footer.n_info_streams; ++i){
		stats_basic[21] += this->info_containers[i];
		stats_info[this->footer.info_offsets[i].data_header.global_key] += this->info_containers[i];
	}

	for(U32 i = 0; i < this->footer.n_format_streams; ++i){
		stats_basic[22] += this->format_containers[i];
		stats_format[this->footer.format_offsets[i].data_header.global_key] += this->format_containers[i];
	}
}

bool VariantBlock::write(std::ostream& stream,
                     import_stats_type& stats_basic,
                     import_stats_type& stats_info,
                     import_stats_type& stats_format)
{
	const U64 begin_pos = stream.tellp();
	this->header.l_offset_footer = this->__determineCompressedSize();
	stream << this->header;
	const U64 start_pos = stream.tellp();
	stats_basic[0].cost_uncompressed += start_pos - begin_pos;

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted){
		this->footer.offset_ppa = this->ppa_manager.header;
		this->footer.offset_ppa.data_header.offset = (U64)stream.tellp() - start_pos;
		stream << this->ppa_manager;
	}

	this->__writeContainer(stream, this->footer.offset_meta_contig,      this->meta_contig_container,    (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_position,    this->meta_positions_container, (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_controllers, this->meta_controller_container,(U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_quality,     this->meta_quality_container,   (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_names,       this->meta_names_container,     (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_refalt,      this->meta_refalt_container,    (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_alleles,     this->meta_alleles_container,   (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_info_id,     this->meta_info_map_ids,        (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_filter_id,   this->meta_filter_map_ids,      (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_meta_format_id,   this->meta_format_map_ids,      (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_helper,        this->gt_support_data_container,(U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_8b,            this->gt_rle8_container,        (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_16b,           this->gt_rle16_container,       (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_32b,           this->gt_rle32_container,       (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_64b,           this->gt_rle64_container,       (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_simple8,       this->gt_simple8_container,     (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_simple16,      this->gt_simple16_container,    (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_simple32,      this->gt_simple32_container,    (U64)stream.tellp() - start_pos);
	this->__writeContainer(stream, this->footer.offset_gt_simple64,      this->gt_simple64_container,    (U64)stream.tellp() - start_pos);

	for(U32 i = 0; i < this->footer.n_info_streams; ++i)
		this->__writeContainer(stream, this->footer.info_offsets[i], this->info_containers[i], (U64)stream.tellp() - start_pos);

	for(U32 i = 0; i < this->footer.n_format_streams; ++i)
		this->__writeContainer(stream, this->footer.format_offsets[i], this->format_containers[i], (U64)stream.tellp() - start_pos);

	// writing footer
	assert(this->header.l_offset_footer == (U64)stream.tellp() - start_pos);

	// Update stats
	this->updateOutputStatistics(stats_basic, stats_info, stats_format);

	return(stream.good());
}

bool VariantBlock::operator+=(meta_entry_type& meta_entry){
	// Meta positions
	this->meta_positions_container.Add((S32)meta_entry.position);
	++this->meta_positions_container;

	// Contig ID
	this->meta_contig_container.Add((S32)meta_entry.contigID);
	++this->meta_contig_container;

	// Ref-alt data
	if(meta_entry.usePackedRefAlt()){ // Is simple SNV and possible extra case when <NON_REF> in gVCF
		meta_entry.controller.alleles_packed = true;
		const BYTE ref_alt = meta_entry.packRefAltByte();
		this->meta_refalt_container.AddLiteral(ref_alt);
		++this->meta_refalt_container;
	}
	// add complex
	else {
		// Special encoding
		for(U32 i = 0; i < meta_entry.n_alleles; ++i){
			// Write out allele
			this->meta_alleles_container.AddLiteral((U16)meta_entry.alleles[i].l_allele);
			this->meta_alleles_container.AddCharacter(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
		}
		++this->meta_alleles_container; // update before to not trigger
		this->meta_alleles_container.addStride(meta_entry.n_alleles);
	}

	// Quality
	this->meta_quality_container.Add(meta_entry.quality);
	++this->meta_quality_container;

	// Variant name
	this->meta_names_container.addStride(meta_entry.name.size());
	this->meta_names_container.AddCharacter(meta_entry.name);
	++this->meta_names_container;

	// Tachyon pattern identifiers
	this->meta_info_map_ids.Add(meta_entry.info_pattern_id);
	this->meta_format_map_ids.Add(meta_entry.format_pattern_id);
	this->meta_filter_map_ids.Add(meta_entry.filter_pattern_id);
	++this->meta_info_map_ids;
	++this->meta_format_map_ids;
	++this->meta_filter_map_ids;

	// Check if all variants are of length 1 (as in all alleles are SNVs)
	bool all_snv = true;
	for(U32 i = 0; i < meta_entry.n_alleles; ++i){
		if(meta_entry.alleles[i].l_allele != 1) all_snv = false;
	}
	meta_entry.controller.all_snv = all_snv;

	// Controller
	this->meta_controller_container.AddLiteral((U16)meta_entry.controller.toValue()); // has been overloaded
	++this->meta_controller_container;

	return true;
}

}
}
