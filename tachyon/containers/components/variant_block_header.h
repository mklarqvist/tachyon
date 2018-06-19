#ifndef INDEX_INDEXBLOCKENTRY_H_
#define INDEX_INDEXBLOCKENTRY_H_

#include <fstream>
#include <bitset>

#include "../hash_container.h"
#include "../../io/basic_buffer.h"
#include "../components/data_container_header_controller.h"
#include "data_block_bitvector.h"
#include "data_container_header.h"

namespace tachyon{
namespace containers{

/** @brief Controller flags for an IndexBlockEntry
 * This structure is for internal use only and describes
 * various internal states as flags.
 */
struct VariantBlockHeaderController{
	typedef VariantBlockHeaderController self_type;

public:
	VariantBlockHeaderController():
		hasGT(0),
		hasGTPermuted(0),
		anyEncrypted(0),
		unused(0)
	{}
	~VariantBlockHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& controller){
		const U16 c = controller.hasGT |
                      controller.hasGTPermuted << 1 |
                      controller.anyEncrypted  << 2 |
                      controller.unused        << 3;

		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(U16));
		return(stream);
	}

public:
	U16 hasGT:         1,  // This block has GT FORMAT data
		hasGTPermuted: 1,  // have the GT fields been permuted
		anyEncrypted:  1,  // any data encrypted
		unused:        13; // reserved for future use
};

/** @brief Fixed-sized components of an IndexBlockEntry
 * For internal use only. This is a subcomponent of fixed
 * size describing:
 * 1) Contig, minimum and maximum genomic coordinates
 * 2) Offsets into the containers
 * 3) Number of containers and ID patterns
 * 4) Controller flags
 */
struct VariantBlockHeader{
private:
	typedef VariantBlockHeader           self_type;
	typedef VariantBlockHeaderController controller_type;
	typedef DataContainerHeader          header_type;

public:
	VariantBlockHeader();
	~VariantBlockHeader();

	inline const U32& size(void) const{ return(this->n_variants); }
	inline const S32& getContigID(void) const{ return(this->contigID); }
	inline const S64& getMinPosition(void) const{ return(this->minPosition); }
	inline const S64& getMaxPosition(void) const{ return(this->maxPosition); }
	inline U64& getBlockID(void){ return(this->blockID); }
	inline const U64& getBlockID(void) const{ return(this->blockID); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry){
		stream.write(reinterpret_cast<const char*>(&entry.l_offset_footer),   sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.blockID),           sizeof(U64));
		stream << entry.controller;
		stream.write(reinterpret_cast<const char*>(&entry.contigID),          sizeof(U32));
		stream.write(reinterpret_cast<const char*>(&entry.minPosition),       sizeof(S64));
		stream.write(reinterpret_cast<const char*>(&entry.maxPosition),       sizeof(S64));
		stream.write(reinterpret_cast<const char*>(&entry.n_variants),        sizeof(U32));

		return(stream);
	}

	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry){
		stream.read(reinterpret_cast<char*>(&entry.l_offset_footer),   sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.blockID),           sizeof(U64));
		stream >> entry.controller;
		stream.read(reinterpret_cast<char*>(&entry.contigID),          sizeof(U32));
		stream.read(reinterpret_cast<char*>(&entry.minPosition),       sizeof(S64));
		stream.read(reinterpret_cast<char*>(&entry.maxPosition),       sizeof(S64));
		stream.read(reinterpret_cast<char*>(&entry.n_variants),        sizeof(U32));

		return(stream);
	}

	void reset(void){
		this->l_offset_footer    = 0;
		this->blockID            = 0;
		this->controller.clear();
		this->contigID           = -1;
		this->minPosition        = 0;
		this->maxPosition        = 0;
		this->n_variants         = 0;
	}

public:
	// allows jumping to the next block when streaming
	// over the file and not using the index
	// EOF marker is at this position - sizeof(EOF marker)
	U32 l_offset_footer;
	U64 blockID;        // block identifier in the form of a random hash
	controller_type controller;
	S32 contigID;       // contig identifier
	S64 minPosition;    // minimum coordinate in this block
	S64 maxPosition;    // maximum coordinate in this block
	U32 n_variants;     // number of variants in this block
};

}
}

#endif /* INDEX_IndexBLOCKENTRY_H_ */
