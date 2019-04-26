#include "basic_reader.h"
#include "utility.h"

namespace tachyon{
namespace io{

BasicReader::BasicReader() :
	filesize_(0),
	block_size_(65536),
	capacity_(this->block_size_*2),
	end_(0),
	buffer_(new type[this->capacity_])
{}

BasicReader::BasicReader(std::string input) :
	filename_(input),
	filesize_(0),
	block_size_(65536),
	capacity_(this->block_size_*2),
	end_(0),
	buffer_(new type[this->capacity_])
{}

BasicReader::BasicReader(std::string input, const size_t block_size) :
	filename_(input),
	filesize_(0),
	block_size_(block_size),
	capacity_(this->block_size_*2),
	end_(0),
	buffer_(new type[this->capacity_])
{}

BasicReader::BasicReader(const self_type& other) :
	filename_(other.filename_),
	filesize_(other.filesize_),
	block_size_(other.block_size_),
	capacity_(other.capacity_),
	end_(other.end_),
	buffer_(new type[this->capacity_])
{
	memcpy(this->buffer_, other.buffer_, other.end_);
	this->open();
}

bool BasicReader::open(std::string filename) {
	// If filename is empty
	if (filename.size() == 0)
		return false;

	// Set filename
	this->filename_ = filename;

	// Open the file
	return(this->open());
}

bool BasicReader::open(void) {
	// Check that filename is set
	if (this->filename_.size() == 0) {
		std::cerr << utility::timestamp("ERROR", "IO") << "No input file given..." << std::endl;
		return false;
	}

	// Open stream at the end of the file
	this->stream_.open(this->filename_, std::ios::binary | std::ios::ate);
	if (!this->good()) {
		std::cerr << utility::timestamp("ERROR", "IO") << "Failed to open file..." << std::endl;
		return false;
	}

	// Set file size
	this->filesize_ = this->stream_.tellg();
	this->stream_.seekg(0);

	// If filesize is 0 return
	if (this->filesize_ <= 0) {
		std::cerr << utility::timestamp("ERROR", "IO") << "File size is 0..." << std::endl;
		return false;
	}

	// Reset buffer pointer to 0
	this->end_ = 0;

	//if (!SILENT)
	//	std::cerr << utility::timestamp("LOG", "IO") << "Opened file: " << this->filename_ << " (" << this->filesize_ << " b)..." << std::endl;

	return(true);
}

void BasicReader::close(void) { this->stream_.close(); }

bool BasicReader::read(void) {
	if (!this->good())
		return false;

	this->end_ = 0;
	this->stream_.read(&this->buffer_[0], this->block_size_);
	this->end_ = this->stream_.gcount();

	if (this->stream_.gcount() > 0)
		return true;
	else
		return false;
}

bool BasicReader::read(const uint32_t length) {
	if (!this->good())
		return false;

	this->end_ = 0;
	this->stream_.read(&this->buffer_[0], length*sizeof(type));
	this->end_ = this->stream_.gcount();

	if (this->stream_.gcount() > 0)
		return true;
	else
		return false;
}

bool BasicReader::readAppend(const uint32_t length) {
	if (!this->good())
		return false;

	this->stream_.read(&this->buffer_[this->end_], length*sizeof(type));
	this->end_ += this->stream_.gcount();

	if (this->stream_.gcount() > 0)
		return true;
	else
		return false;
}

bool BasicReader::getLine(void) { // Read until finding a new line into buffer
	if (!this->good())
		return false;

	this->stream_.getline(&this->buffer_[this->end_], this->capacity_ - this->end_);
	this->end_ += this->stream_.gcount();

	if (this->stream_.eof()) {
		//std::cerr << Tachyon::utility::timestamp("LOG", "IO") << "EOF found check bit: " << this->stream_.eof() << std::endl;
		return false;
	}

	if ((this->stream_.fail()) && (this->capacity_ - this->end_ - 1 == 0)) {
		if (!SILENT)
			std::cerr << utility::timestamp("LOG", "IO") << "Resizing buffer from: " << this->capacity_ << " -> " << this->capacity_*2 << std::endl;

		this->resize();
		this->stream_.clear();
		this->getLine();
	}

	if (this->stream_.fail()) {
		std::cerr << utility::timestamp("ERROR", "IO") << "IO stream failed!" << std::endl;
		std::cerr << utility::timestamp("ERROR", "IO") << "Stream position: " << this->stream_.tellg() << std::endl;
		return false;
	}

	return true;
}

bool BasicReader::getLine(std::string& data) { // Read until finding a new line into string
	if (!this->good())
		return false;

	std::getline(this->stream_, data);
	if (this->stream_.fail())
		return false;

	else return true;
}

}
}
