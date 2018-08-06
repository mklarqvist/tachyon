#include "variant_block.h"
#include "algorithm/compression/compression_container.h"
#include "support/helpers.h"
#include "support/type_definitions.h"

namespace tachyon{
namespace containers{

VariantBlock::VariantBlock() :
	n_info_c_allocated(0),
	n_format_c_allocated(0),
	base_containers(new container_type[YON_BLK_N_STATIC]),
	info_containers(new container_type[0]),
	format_containers(new container_type[0]),
	gt_ppa(nullptr),
	end_block_(0),
	start_compressed_data_(0),
	end_compressed_data_(0)
{
	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);
	this->footer_support.resize(65536);
}

VariantBlock::VariantBlock(const uint16_t n_info, const uint16_t n_format) :
	n_info_c_allocated(n_info),
	n_format_c_allocated(n_format),
	base_containers(new container_type[YON_BLK_N_STATIC]),
	info_containers(new container_type[n_info]),
	format_containers(new container_type[n_format]),
	gt_ppa(nullptr),
	end_block_(0),
	start_compressed_data_(0),
	end_compressed_data_(0)
{
	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);
	this->footer_support.resize(65536);
}

VariantBlock::~VariantBlock(){
	delete [] this->base_containers;
	delete [] this->info_containers;
	delete [] this->format_containers;
	delete this->gt_ppa;
}

void VariantBlock::clear(void){
	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i)
		this->base_containers[i].reset();

	for(U32 i = 0; i < this->footer.n_info_streams; ++i)
		this->info_containers[i].reset();

	for(U32 i = 0; i < this->footer.n_format_streams; ++i)
		this->format_containers[i].reset();

	this->base_containers[YON_BLK_ALLELES].SetType(YON_TYPE_STRUCT);
	this->base_containers[YON_BLK_CONTROLLER].SetType(YON_TYPE_16B);
	this->base_containers[YON_BLK_REFALT].SetType(YON_TYPE_8B);

	this->end_block_             = 0;
	this->start_compressed_data_ = 0;
	this->end_compressed_data_   = 0;

	this->header.reset();
	this->footer.reset();
	this->footer_support.reset();

	if(this->gt_ppa != nullptr)
		this->gt_ppa->reset();
}

void VariantBlock::resize(const U32 s){
	if(s == 0) return;

	for(U32 i = 0; i < YON_BLK_N_STATIC; ++i)
		this->base_containers[i].resize(s);

	for(U32 i = 0; i < n_info_c_allocated; ++i){
		this->info_containers[i].resize(s);
	}

	for(U32 i = 0; i < n_format_c_allocated; ++i){
		this->format_containers[i].resize(s);
	}
}

void VariantBlock::UpdateContainers(void){
	this->base_containers[YON_BLK_CONTIG].UpdateContainer();
	this->base_containers[YON_BLK_POSITION].UpdateContainer();
	this->base_containers[YON_BLK_REFALT].UpdateContainer(false, false);
	this->base_containers[YON_BLK_QUALITY].UpdateContainer();
	this->base_containers[YON_BLK_NAMES].UpdateContainer();
	this->base_containers[YON_BLK_ALLELES].UpdateContainer(false, false);
	this->base_containers[YON_BLK_ID_FILTER].UpdateContainer();
	this->base_containers[YON_BLK_ID_FORMAT].UpdateContainer();
	this->base_containers[YON_BLK_ID_INFO].UpdateContainer();
	this->base_containers[YON_BLK_GT_SUPPORT].UpdateContainer();
	this->base_containers[YON_BLK_GT_PLOIDY].UpdateContainer();

	this->base_containers[YON_BLK_CONTROLLER].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_INT8].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT16].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT32].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_INT64].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT8].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT16].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT32].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_S_INT64].UpdateContainer(false, true);
	this->base_containers[YON_BLK_GT_N_INT8].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT16].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT32].UpdateContainer(false, false);
	this->base_containers[YON_BLK_GT_N_INT64].UpdateContainer(false, false);

	for(U32 i = 0; i < this->footer.n_info_streams; ++i){
		assert(this->info_containers[i].header.data_header.stride != 0);
		this->info_containers[i].UpdateContainer();
	}

	for(U32 i = 0; i < this->footer.n_format_streams; ++i){
		assert(this->format_containers[i].header.data_header.stride != 0);
		this->format_containers[i].UpdateContainer();
	}
}

bool VariantBlock::ReadHeaderFooter(std::ifstream& stream){
	if(!stream.good()){
		std::cerr << utility::timestamp("ERROR") << "File stream is corrupted..." << std::endl;
		return false;
	}

	stream >> this->header; // load header
	this->start_compressed_data_ = (U64)stream.tellg(); // start of compressed data
	stream.seekg(this->start_compressed_data_ + this->header.l_offset_footer); // seek to start of footer
	this->end_compressed_data_   = stream.tellg(); // end of compressed data

	assert(stream.good());

	U32 footer_uLength = 0;
	U32 footer_cLength = 0;
	uint8_t footer_crc[MD5_DIGEST_LENGTH];
	utility::DeserializePrimitive(footer_uLength, stream);
	utility::DeserializePrimitive(footer_cLength, stream);
	stream.read(reinterpret_cast<char*>(&footer_crc[0]), MD5_DIGEST_LENGTH);

	this->footer_support.resize(footer_cLength);
	stream.read(this->footer_support.buffer_data.data(), footer_cLength);
	this->footer_support.buffer_data.n_chars = footer_cLength;
	this->footer_support.buffer_data_uncompressed.resize(footer_uLength);
	this->footer_support.buffer_data_uncompressed.n_chars      = footer_uLength;
	this->footer_support.header.data_header.controller.encoder = YON_ENCODE_ZSTD;
	this->footer_support.header.data_header.cLength            = footer_cLength;
	this->footer_support.header.data_header.uLength            = footer_uLength;
	memcpy(&this->footer_support.header.data_header.crc[0], &footer_crc[0], MD5_DIGEST_LENGTH);

	// Assert end-of-block marker
	U64 eof_marker;
	utility::DeserializePrimitive(eof_marker, stream);
	assert(eof_marker == constants::TACHYON_BLOCK_EOF);
	this->end_block_ = stream.tellg(); // end-of-block offset
	stream.seekg(this->start_compressed_data_);
	return(stream.good());
}

bool VariantBlock::read(std::ifstream& stream){
	if(this->header.controller.hasGTPermuted && this->header.controller.hasGT){
		stream.seekg(this->start_compressed_data_ + this->footer.offsets[YON_BLK_PPA].data_header.offset);
		stream >> this->base_containers[YON_BLK_PPA];
	}

	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i)
		this->LoadContainer(stream, this->footer.offsets[i], this->base_containers[i]);

	// Load all INFO
	delete [] this->info_containers;
	this->info_containers = new container_type[this->footer.n_info_streams];
	this->n_info_c_allocated = this->footer.n_info_streams;
	if(this->footer.n_info_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.info_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_info_streams; ++i)
			this->LoadContainer(stream, this->footer.info_offsets[i], this->info_containers[i]);

	}

	// Load all FORMAT
	delete [] this->format_containers;
	this->format_containers = new container_type[this->footer.n_format_streams];
	this->n_format_c_allocated = this->footer.n_format_streams;
	if(this->footer.n_format_streams){
		stream.seekg(this->start_compressed_data_ + this->footer.format_offsets[0].data_header.offset);
		for(U32 i = 0; i < this->footer.n_format_streams; ++i)
			this->LoadContainer(stream, this->footer.format_offsets[i], this->format_containers[i]);

		// EOF assertion
		assert(this->end_compressed_data_ == (U64)stream.tellg());
	}

	stream.seekg(this->end_block_); // seek to end-of-block
	return(true);
}

U64 VariantBlock::DetermineCompressedSize(void) const{
	U64 total = 0;
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		total += this->base_containers[YON_BLK_PPA].GetObjectSize();

	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i)              total += this->base_containers[i].GetObjectSize();
	for(U32 i = 0; i < this->footer.n_info_streams; ++i)   total += this->info_containers[i].GetObjectSize();
	for(U32 i = 0; i < this->footer.n_format_streams; ++i) total += this->format_containers[i].GetObjectSize();

	return(total);
}

void VariantBlock::UpdateOutputStatistics(import_stats_type& stats_basic, import_stats_type& stats_info, import_stats_type& stats_format){
	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		stats_basic[1] += this->base_containers[YON_BLK_PPA];

	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i)
		stats_basic[i+1] += this->base_containers[i];

	for(U32 i = 0; i < this->footer.n_info_streams; ++i){
		stats_basic[22] += this->info_containers[i];
		stats_info[this->footer.info_offsets[i].data_header.global_key] += this->info_containers[i];
	}

	for(U32 i = 0; i < this->footer.n_format_streams; ++i){
		stats_basic[23] += this->format_containers[i];
		stats_format[this->footer.format_offsets[i].data_header.global_key] += this->format_containers[i];
	}
}

bool VariantBlock::write(std::ostream& stream,
                     import_stats_type& stats_basic,
                     import_stats_type& stats_info,
                     import_stats_type& stats_format)
{
	if(stream.good() == false){
		return false;
	}

	const U64 begin_pos = stream.tellp();
	this->header.l_offset_footer = this->DetermineCompressedSize();
	stream << this->header;
	const U64 start_pos = stream.tellp();
	stats_basic[0].cost_uncompressed += start_pos - begin_pos;

	if(this->header.controller.hasGT && this->header.controller.hasGTPermuted)
		this->WriteContainer(stream, this->footer.offsets[YON_BLK_PPA], this->base_containers[YON_BLK_PPA], (U64)stream.tellp() - start_pos);

	// Start at offset 1 because offset 0 is encoding for the genotype
	// permutation array that is handled differently.
	for(U32 i = 1; i < YON_BLK_N_STATIC; ++i)
		this->WriteContainer(stream, this->footer.offsets[i], this->base_containers[i], (U64)stream.tellp() - start_pos);

	for(U32 i = 0; i < this->footer.n_info_streams; ++i)
		this->WriteContainer(stream, this->footer.info_offsets[i], this->info_containers[i], (U64)stream.tellp() - start_pos);

	for(U32 i = 0; i < this->footer.n_format_streams; ++i)
		this->WriteContainer(stream, this->footer.format_offsets[i], this->format_containers[i], (U64)stream.tellp() - start_pos);

	// Assert that the written amount equals the expected amount.
	assert(this->header.l_offset_footer == (U64)stream.tellp() - start_pos);

	// Update stats
	this->UpdateOutputStatistics(stats_basic, stats_info, stats_format);

	return(stream.good());
}

bool VariantBlock::operator+=(meta_entry_type& meta_entry){
	// Meta positions
	this->base_containers[YON_BLK_POSITION].Add((S32)meta_entry.position);
	++this->base_containers[YON_BLK_POSITION];

	// Contig ID
	this->base_containers[YON_BLK_CONTIG].Add((S32)meta_entry.contigID);
	++this->base_containers[YON_BLK_CONTIG];

	// Ref-alt data
	if(meta_entry.UsePackedRefAlt()){ // Is simple SNV and possible extra case when <NON_REF> in gVCF
		meta_entry.controller.alleles_packed = true;
		const BYTE ref_alt = meta_entry.PackRefAltByte();
		this->base_containers[YON_BLK_REFALT].AddLiteral(ref_alt);
		++this->base_containers[YON_BLK_REFALT];
	}
	// add complex
	else {
		// Special encoding
		for(U32 i = 0; i < meta_entry.n_alleles; ++i){
			// Write out allele
			this->base_containers[YON_BLK_ALLELES].AddLiteral((U16)meta_entry.alleles[i].l_allele);
			this->base_containers[YON_BLK_ALLELES].AddCharacter(meta_entry.alleles[i].allele, meta_entry.alleles[i].l_allele);
		}
		++this->base_containers[YON_BLK_ALLELES]; // update before to not trigger
		this->base_containers[YON_BLK_ALLELES].AddStride(meta_entry.n_alleles);
	}

	// Quality
	this->base_containers[YON_BLK_QUALITY].Add(meta_entry.quality);
	++this->base_containers[YON_BLK_QUALITY];

	// Variant name
	this->base_containers[YON_BLK_NAMES].AddStride(meta_entry.name.size());
	this->base_containers[YON_BLK_NAMES].AddCharacter(meta_entry.name);
	++this->base_containers[YON_BLK_NAMES];

	// Tachyon pattern identifiers
	this->base_containers[YON_BLK_ID_INFO].Add(meta_entry.info_pattern_id);
	this->base_containers[YON_BLK_ID_FORMAT].Add(meta_entry.format_pattern_id);
	this->base_containers[YON_BLK_ID_FILTER].Add(meta_entry.filter_pattern_id);
	++this->base_containers[YON_BLK_ID_INFO];
	++this->base_containers[YON_BLK_ID_FORMAT];
	++this->base_containers[YON_BLK_ID_FILTER];

	// Check if all variants are of length 1 (as in all alleles are SNVs)
	bool all_snv = true;
	for(U32 i = 0; i < meta_entry.n_alleles; ++i){
		if(meta_entry.alleles[i].size() != 1) all_snv = false;
	}
	meta_entry.controller.all_snv = all_snv;

	// Controller
	this->base_containers[YON_BLK_CONTROLLER].AddLiteral((U16)meta_entry.controller.toValue()); // has been overloaded
	++this->base_containers[YON_BLK_CONTROLLER];

	// Ploidy
	this->base_containers[YON_BLK_GT_PLOIDY].Add(meta_entry.n_base_ploidy);
	++this->base_containers[YON_BLK_GT_PLOIDY];

	return true;
}

}
}
