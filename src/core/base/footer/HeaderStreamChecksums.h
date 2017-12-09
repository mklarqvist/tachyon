#ifndef CORE_BASE_FOOTER_FOOTERSTREAMCHECKSUMS_H_
#define CORE_BASE_FOOTER_FOOTERSTREAMCHECKSUMS_H_

namespace Tachyon{
namespace Core{

struct FooterStreamChecksums{
private:
	typedef FooterStreamChecksums self_type;

public:
	FooterStreamChecksums(void) : n_width(0), uncompressed_checksum(nullptr), compressed_checksum(nullptr){}
	~FooterStreamChecksums(){
		delete [] this->uncompressed_checksum;
		delete [] this->compressed_checksum;
	}

	bool setHashes(const U16& n_width, const BYTE* uncompressed, const BYTE* compressed){
		if(!this->setUncompressedHash(n_width, uncompressed)){
			return false;
		}
		if(!this->setCompressedHash(n_width, compressed)){
			return false;
		}
		return true;
	}

	bool setUncompressedHash(const U16& n_width, const BYTE* uncompressed){
		if(this->n_width != n_width && this->n_width != 0){
			std::cerr << "cannot provide two different sum lengths" << std::endl;
			return false;
		}
		delete [] this->uncompressed_checksum;
		this->uncompressed_checksum = new BYTE[n_width];
		this->n_width = n_width;
		memcpy(this->uncompressed_checksum, uncompressed, n_width);
		return true;
	}

	bool setCompressedHash(const U16& n_width, const BYTE* compressed){
		if(this->n_width != n_width && this->n_width != 0){
			std::cerr << "cannot provide two different sum lengths" << std::endl;
			return false;
		}
		delete [] this->compressed_checksum;
		this->compressed_checksum = new BYTE[n_width];
		this->n_width = n_width;
		memcpy(this->compressed_checksum, compressed, n_width);
		return true;
	}

private:
	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.n_width), sizeof(U16));
		stream.write((char*)&entry.uncompressed_checksum[0], entry.n_width);
		stream.write((char*)&entry.compressed_checksum[0],   entry.n_width);
		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.n_width), sizeof(U16));
		entry.uncompressed_checksum = new BYTE[entry.n_width];
		stream.read((char*)&entry.uncompressed_checksum[0], entry.n_width);
		entry.compressed_checksum = new BYTE[entry.n_width];
		stream.read((char*)&entry.compressed_checksum[0], entry.n_width);
		return(stream);
	}

public:
	U16 n_width;
	BYTE* uncompressed_checksum;
	BYTE*   compressed_checksum;
};

}
}

#endif /* CORE_BASE_FOOTER_FOOTERSTREAMCHECKSUMS_H_ */
