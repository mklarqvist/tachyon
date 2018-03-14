#ifndef IO_BAM_BAMREADER_H_
#define IO_BAM_BAMREADER_H_

#include <cassert>

#include "BCFEntry.h"
#include "../compression/BGZFController.h"


namespace tachyon {
namespace bam {

class BAMReader{
private:
    typedef BAMReader          self_type;
    typedef io::BasicBuffer    buffer_type;
    typedef io::BGZFController bgzf_controller_type;

public:
    BAMReader();
    BAMReader(const std::string& file_name);
	~BAMReader();

	bool open(const std::string input){
		this->file_name = input;
		return(this->open());
	}

	bool open(void){
		if(this->file_name.size() == 0)
			return false;

		this->stream.open(this->file_name, std::ios::binary | std::ios::in | std::ios::ate);
		if(!this->stream.good()){
			std::cerr << utility::timestamp("ERROR", "BAM") << "Failed to open file: " << this->file_name << std::endl;
			return false;
		}

		this->filesize = this->stream.tellg();
		this->stream.seekg(0);

		if(!this->stream.good()){
			std::cerr << utility::timestamp("ERROR", "BAM") << "Bad stream!" << std::endl;
			return false;
		}

		if(!this->nextBlock()){
			std::cerr << utility::timestamp("ERROR","BAM") << "Failed to get first block!" << std::endl;
			return false;
		}

		if(!this->parseHeader()){
			std::cerr << utility::timestamp("ERROR","BAM") << "Failed to parse header!" << std::endl;
			return false;
		}

		return true;
	}

	bool nextBlock(void){
		// Stream died
		if(!this->stream.good()){
			std::cerr << utility::timestamp("ERROR", "BAM") << "Stream died!" << std::endl;
			return false;
		}

		// EOF
		if(this->stream.tellg() == this->filesize){
			return false;
		}

		if(!this->bgzf_controller.InflateBlock(this->stream, this->buffer)){
			return false;
		}

		// Reset buffer
		this->buffer.reset();
		this->current_pointer = 0;
		this->b_data_read += this->bgzf_controller.buffer.size();

		return true;
	}

private:
	bool parseHeader(void){
		if(this->bgzf_controller.buffer.size() == 0){
			std::cerr << utility::timestamp("ERROR","BAM") << "No buffer!" << std::endl;
			return false;
		}

		if(this->bgzf_controller.buffer.size() < 5){
			std::cerr << utility::timestamp("ERROR","BAM") << "Corrupted header!" << std::endl;
			return false;
		}

		if(strncmp(&this->bgzf_controller.buffer.buffer[0], "BAM\2\2", 5) != 0){
			std::cerr << utility::timestamp("ERROR","BAM") << "Failed to validate MAGIC" << std::endl;
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

		return true;
	}

public:
	std::string          file_name;
	std::ifstream        stream;
	U64                  filesize;
	U32                  current_pointer;
	S32                  map_gt_id;
	buffer_type          buffer;
	buffer_type          header_buffer;
	bgzf_controller_type bgzf_controller;
	U64                  b_data_read;
};

}
}



#endif /* IO_BAM_BAMREADER_H_ */
