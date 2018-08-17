#ifndef INDEX_INDEXBLOCKENTRY_H_
#define INDEX_INDEXBLOCKENTRY_H_

#include <fstream>
#include <bitset>

#include "data_container_header.h"
#include "io/basic_buffer.h"
#include "data_container_header_controller.h"

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
	inline const S32& GetContigID(void) const{ return(this->contigID); }
	inline const S64& GetMinPosition(void) const{ return(this->minPosition); }
	inline const S64& GetMaxPosition(void) const{ return(this->maxPosition); }
	inline U64& GetBlockHash(void){ return(this->block_hash); }
	inline const U64& GetBlockHash(void) const{ return(this->block_hash); }

	friend std::ostream& operator<<(std::ostream& stream, const self_type& entry);
	friend std::ifstream& operator>>(std::ifstream& stream, self_type& entry);

	void reset(void);

public:
	// Allows jumping to the next block when streaming.
	// EOF marker is at this position - sizeof(EOF marker)
	U32 l_offset_footer;
	U64 block_hash;     // block identifier in the form of a random hash
	controller_type controller;
	S32 contigID;       // contig identifier
	S64 minPosition;    // minimum coordinate in this block
	S64 maxPosition;    // maximum coordinate in this block
	U32 n_variants;     // number of variants in this block
};

}
}

#endif /* INDEX_IndexBLOCKENTRY_H_ */
