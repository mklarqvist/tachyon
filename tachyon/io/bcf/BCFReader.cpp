#include "BCFReader.h"

namespace tachyon{
namespace bcf{

BCFReader::BCFReader() :
		filesize(0),
		current_pointer(0),
		state(bcf_reader_state::BCF_INIT),
		n_entries(0),
		n_capacity(0),
		n_carry_over(0),
		entries(nullptr)
{}

BCFReader::~BCFReader(){
	if(this->entries != nullptr){
		for(U32 i = 0; i < this->n_capacity; ++i)
			delete this->entries[i];

		delete this->entries;
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

	return true;
}

bool BCFReader::nextVariant(BCFEntry& entry){
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

	if(!entry.parse()){
		std::cerr << "parse error" << std::endl;
		exit(1);
	}

	return true;
}

bool BCFReader::getVariants(const U32 n_variants, const double bp_window, bool across_contigs){
	if(this->n_entries == 0 && this->n_carry_over == 0)
		this->entries = new entry_type*[n_variants + 1];


	for(U32 i = 0; i < this->n_entries; ++i)
		delete this->entries[i];

	// If there is any carry over
	U32 firstPos = 0;
	S32 firstContig = -1;
	if(this->n_carry_over == 1){
		this->entries[0] = this->entries[this->n_entries];
		this->n_entries = 1;
		firstPos = this->entries[0]->body->POS;
		firstContig = this->entries[0]->body->CHROM;
		//std::cerr << utility::timestamp("LOG", "SWITCHING") << firstPos << '\t' << firstContig << std::endl;
	} else {
		// Only set this to 0 if there is no carry
		// over data from the previous cycle
		this->n_entries = 0;
	}

	// Entries
	//this->entries = new entry_type*[n_variants];
	this->n_carry_over = 0;

	// EOF
	if(this->state == bcf_reader_state::BCF_EOF)
		return false;


	for(U32 i = 0; i < n_variants; ++i){
		if(this->current_pointer == this->bgzf_controller.buffer.size()){
			if(!this->nextBlock()){
				return(this->size() > 0);
			}
		}

		this->entries[this->n_entries] = new entry_type;

		if(this->current_pointer + 8 > this->bgzf_controller.buffer.size()){
			const S32 partial = (S32)this->bgzf_controller.buffer.size() - this->current_pointer;
			this->entries[this->n_entries]->add(&this->bgzf_controller.buffer[this->current_pointer], this->bgzf_controller.buffer.size() - this->current_pointer);
			if(!this->nextBlock()){
				std::cerr << utility::timestamp("ERROR","BCF") << "Failed to get next block in partial" << std::endl;
				return false;
			}

			this->entries[this->n_entries]->add(&this->bgzf_controller.buffer[0], 8 - partial);
			this->current_pointer = 8 - partial;
		} else {
			this->entries[this->n_entries]->add(&this->bgzf_controller.buffer[this->current_pointer], 8);
			this->current_pointer += 8;
		}

		U64 remainder = this->entries[this->n_entries]->sizeBody();
		if(remainder > this->entries[this->n_entries]->capacity())
			this->entries[this->n_entries]->resize(remainder + 1024);


		while(remainder > 0){
			if(this->current_pointer + remainder > this->bgzf_controller.buffer.size()){
				this->entries[this->n_entries]->add(&this->bgzf_controller.buffer[this->current_pointer], this->bgzf_controller.buffer.size() - this->current_pointer);
				remainder -= this->bgzf_controller.buffer.size() - this->current_pointer;
				if(!this->nextBlock()){
					std::cerr << utility::timestamp("ERROR","BCF") << "Failed to get next block in partial" << std::endl;
					return false;
				}

			} else {
				this->entries[this->n_entries]->add(&this->bgzf_controller.buffer[this->current_pointer], remainder);
				this->current_pointer += remainder;
				remainder = 0;
				break;
			}
		}

		// Interpret char stream
		this->entries[this->n_entries]->parse();
		//std::cerr << this->entries[this->n_entries]->body->CHROM << ':' << this->entries[this->n_entries]->body->POS+1 << std::endl;


		// Check position
		if(this->n_entries == 0){
			firstPos = this->entries[0]->body->POS;
			firstContig = this->entries[0]->body->CHROM;
		}


		// Make sure that data does not span over
		// multiple CHROM fields
		// Note: this should be toggleable as a
		// passable parameter
		// Note: This property is maintainable only
		// when the input file is sorted
		if(!across_contigs){
			if(this->entries[this->n_entries]->body->CHROM != firstContig){
				std::cerr << utility::timestamp("LOG","CONTIG") << "Switch in CHROM: " << firstContig << "->" << this->entries[this->n_entries]->body->CHROM << std::endl;
				//std::cerr << "Last is now: " << this->last().body->CHROM << ':' << this->last().body->POS << std::endl;
				this->n_carry_over = 1;
				return(this->size() > 0);
			}
		}

		// Check break condition for window
		if(this->entries[this->n_entries]->body->POS - firstPos > bp_window){
			//std::cerr << utility::timestamp("LOG","LD") << "Breaking at " << this->n_entries + 1 << " (" << (this->entries[this->n_entries]->body->POS) - firstPos << ")" << std::endl;
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

bool BCFReader::open(const std::string input){
	this->stream.open(input, std::ios::binary | std::ios::in | std::ios::ate);
	if(!this->stream.good()){
		std::cerr << utility::timestamp("ERROR", "BCF") << "Failed to open file: " << input << std::endl;
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

}
}
