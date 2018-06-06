#include "BCFReader.h"

namespace tachyon{
namespace bcf{

BCFReader::BCFReader() :
		filesize(0),
		current_pointer(0),
		map_gt_id(-1),
		state(bcf_reader_state::BCF_INIT),
		n_entries(0),
		n_capacity(0),
		n_carry_over(0),
		entries(nullptr),
		b_data_read(0)
{}

BCFReader::BCFReader(const std::string& file_name) :
		file_name(file_name),
		filesize(0),
		current_pointer(0),
		map_gt_id(-1),
		state(bcf_reader_state::BCF_INIT),
		n_entries(0),
		n_capacity(0),
		n_carry_over(0),
		entries(nullptr),
		b_data_read(0)
{}

BCFReader::~BCFReader(){
	if(this->entries != nullptr){
		for(std::size_t i = 0; i < this->n_entries; ++i)
			((this->entries + i)->~BCFEntry());

		::operator delete[](static_cast<void*>(this->entries));
	}
}


bool BCFReader::nextBlock(void){
	// Stream died
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "BCF") << "Stream died!" << std::endl;
		this->state = bcf_reader_state::BCF_STREAM_ERROR;
		return false;
	}

	// EOF
	if(this->stream.tellg() == this->filesize){
		this->state = bcf_reader_state::BCF_EOF;
		return false;
	}

	if(!this->bgzf_controller.InflateBlock(this->stream, this->buffer)){
		if(this->bgzf_controller.buffer.size() == 0) this->state = bcf_reader_state::BCF_EOF;
		else this->state = bcf_reader_state::BCF_ERROR;
		return false;
	}

	// Reset buffer
	this->buffer.reset();
	this->current_pointer = 0;
	this->state = bcf_reader_state::BCF_OK;
	this->b_data_read += this->bgzf_controller.buffer.size();

	return true;
}

bool BCFReader::nextVariant(reference entry){
	if(this->current_pointer == this->bgzf_controller.buffer.size()){
		if(!this->nextBlock())
			return false;
	}

	if(this->current_pointer + 8 > this->bgzf_controller.buffer.size()){
		const S32 partial = (S32)this->bgzf_controller.buffer.size() - this->current_pointer;
		entry.add(&this->bgzf_controller.buffer[this->current_pointer], this->bgzf_controller.buffer.size() - this->current_pointer);
		if(!this->nextBlock()){
			std::cerr << utility::timestamp("ERROR","BCF") << "Failed to get next block in partial" << std::endl;
			return false;
		}

		entry.add(&this->bgzf_controller.buffer[0], 8 - partial);
		this->current_pointer = 8 - partial;
	} else {
		entry.add(&this->bgzf_controller.buffer[this->current_pointer], 8);
		this->current_pointer += 8;
	}

	U64 remainder = entry.sizeBody();
	while(remainder){
		if(this->current_pointer + remainder > this->bgzf_controller.buffer.size()){
			entry.add(&this->bgzf_controller.buffer[this->current_pointer], this->bgzf_controller.buffer.size() - this->current_pointer);
			remainder -= this->bgzf_controller.buffer.size() - this->current_pointer;
			if(!this->nextBlock())
				return false;

		} else {
			entry.add(&this->bgzf_controller.buffer[this->current_pointer], remainder);
			this->current_pointer += remainder;
			remainder = 0;
			break;
		}
	}

	if(!entry.parse(this->header.samples)){
		std::cerr << "parse error" << std::endl;
		exit(1);
	}

	if(this->entries[this->n_entries].body->n_fmt > 0 && this->map_gt_id != -1){
		if(this->entries[this->n_entries].formatID[0].mapID == this->map_gt_id){
			this->entries[this->n_entries].hasGenotypes = true;
			entry.assessGenotypes(this->header.samples);
		}
	}

	return true;
}

bool BCFReader::getVariants(const U32 n_variants, const double bp_window, bool across_contigs){
	S64 firstPos    = 0;
	S32 firstContig = -1;
	// If there is any carry over
	if(this->n_carry_over == 1){
		value_type last(this->entries[this->n_entries]);
		if(last.body->n_fmt > 0 && this->map_gt_id != -1){
			if(last.formatID[0].mapID == this->map_gt_id)
				last.hasGenotypes = true;
		}

		for(std::size_t i = 0; i <= this->n_entries; ++i)
			((this->entries + i)->~BCFEntry());

		if(n_variants + 1 > this->capacity()){
			::operator delete[](static_cast<void*>(this->entries));
			this->entries    = static_cast<pointer>(::operator new[]((n_variants + 1)*sizeof(value_type)));
			this->n_capacity = n_variants + 1;
		}

		new( &this->entries[0] ) value_type( last );
		this->n_entries  = 1;
		firstPos         = this->entries[0].body->POS;
		firstContig      = this->entries[0].body->CHROM;
	}
	// Nothing carried over
	else {
		// Only set this to 0 if there is no carry
		// over data from the previous cycle
		for(std::size_t i = 0; i < this->n_entries; ++i)
			((this->entries + i)->~BCFEntry());

		if(n_variants + 1 > this->capacity()){
			::operator delete[](static_cast<void*>(this->entries));
			this->entries    = static_cast<pointer>(::operator new[]((n_variants + 1)*sizeof(value_type)));
			this->n_capacity = n_variants + 1;
		}

		//delete this->entries;
		this->n_entries  = 0;
	}

	// Entries
	this->n_carry_over = 0;

	// EOF
	if(this->state == bcf_reader_state::BCF_EOF)
		return false;

	U32 retrieved_variants = 0;
	bool is_new = true;

	while(retrieved_variants < n_variants){
	//for(U32 i = 0; i < n_variants; ++i){
		if(this->current_pointer == this->bgzf_controller.buffer.size()){
			if(!this->nextBlock()){
				return(this->size() > 0);
			}
		}

		if(is_new) new( &this->entries[this->n_entries] ) value_type( this->header.samples * 2 );
		if(!this->nextVariant(this->entries[this->n_entries])){
			std::cerr << "failed to get next" << std::endl;
			return false;
		}


		if(this->entries[this->n_entries].gt_support.invariant == true){
			//std::cerr << "not getting " << std::endl;
			//((this->entries + this->n_entries)->~BCFEntry());
			//&this->entries[this->n_entries] = static_cast<pointer>(::operator new[](sizeof(value_type)));
			this->entries[this->n_entries].reset();
			is_new = false;
			continue;
		}
		is_new = true;
		++retrieved_variants;

		// Check position
		if(this->n_entries == 0){
			firstPos    = this->entries[0].body->POS;
			firstContig = this->entries[0].body->CHROM;
		}

		// Make sure that data does not span over
		// multiple CHROM fields
		// Note: This property is maintainable only
		// when the input file is sorted
		if(!across_contigs){
			if(this->entries[this->n_entries].body->CHROM != firstContig){
				//std::cerr << utility::timestamp("LOG","CONTIG") << "Switch in CHROM: " << firstContig << "->" << this->entries[this->n_entries].body->CHROM << std::endl;
				this->n_carry_over = 1;
				return(this->size() > 0);
			}
		}

		// Check break condition for window
		if(this->entries[this->n_entries].body->POS - firstPos > bp_window){
			++this->n_entries;
			break;
		}

		// Increment entries in return block
		++this->n_entries;
	}

	return(this->size() > 0);
}

bool BCFReader::parseHeader(void){
	if(this->bgzf_controller.buffer.size() == 0){
		std::cerr << utility::timestamp("ERROR","BCF") << "No buffer!" << std::endl;
		return false;
	}

	if(this->bgzf_controller.buffer.size() < 5){
		std::cerr << utility::timestamp("ERROR","BCF") << "Corrupted header!" << std::endl;
		return false;
	}

	if(strncmp(&this->bgzf_controller.buffer.buffer[0], "BCF\2\2", 5) != 0){
		std::cerr << utility::timestamp("ERROR","BCF") << "Failed to validate MAGIC" << std::endl;
		return false;
	}

	const U32 l_text = *reinterpret_cast<const U32* const>(&this->bgzf_controller.buffer[5]) + 4;
	this->header_buffer.resize(l_text + 1);

	if(l_text - 5 < this->bgzf_controller.buffer.size()){
		this->header_buffer.Add(&this->bgzf_controller.buffer[5], l_text);
		this->current_pointer = l_text + 5;
	} else {
		U32 head_read = this->bgzf_controller.buffer.size() - 5;
		this->header_buffer.Add(&this->bgzf_controller.buffer[5], this->bgzf_controller.buffer.size() - 5);

		//U32 p = 0;
		while(this->nextBlock()){
			if(head_read + this->bgzf_controller.buffer.size() >= l_text){
				this->header_buffer.Add(&this->bgzf_controller.buffer[0], l_text - head_read);
				this->current_pointer = l_text - head_read;
				break;
			}
			head_read += this->bgzf_controller.buffer.size();
			this->header_buffer.Add(&this->bgzf_controller.buffer[0], this->bgzf_controller.buffer.size());
		}
	}

	if(!this->header.parse(&this->header_buffer[0], this->header_buffer.size())){
		std::cerr << utility::timestamp("ERROR","BCF") << "Failed to parse header!" << std::endl;
		return false;
	}

	return true;
}

bool BCFReader::open(){
	if(this->file_name.size() == 0)
		return false;

	this->stream.open(this->file_name, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "BCF") << "Failed to open file: " << this->file_name << std::endl;
		return false;
	}

	this->filesize = this->stream.tellg();
	this->stream.seekg(0);

	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "BCF") << "Bad stream!" << std::endl;
		return false;
	}

	if(!this->nextBlock()){
		std::cerr << utility::timestamp("ERROR","BCF") << "Failed to get first block!" << std::endl;
		return false;
	}

	if(!this->parseHeader()){
		std::cerr << utility::timestamp("ERROR","BCF") << "Failed to parse header!" << std::endl;
		return false;
	}

	return true;
}

bool BCFReader::open(const std::string input){
	this->file_name = input;
	return(this->open());
}

}
}
