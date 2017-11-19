#ifndef CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_
#define CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_

#include <fstream>

namespace Tachyon{
namespace Core{

// Controller type for stream container
struct StreamContainerHeaderController{
	typedef StreamContainerHeaderController self_type;

public:
	StreamContainerHeaderController() :
		signedness(0),
		mixedStride(0),
		type(0),
		encoder(0),
		uniform(0),
		unused(0)
	{}
	~StreamContainerHeaderController(){}

	inline void clear(){ memset(this, 0, sizeof(U16)); }

	friend std::ofstream& operator<<(std::ofstream& stream, const self_type& controller){
		const U16* c = reinterpret_cast<const U16* const>(&controller);
		stream.write(reinterpret_cast<const char*>(&c), sizeof(U16));
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& controller){
		U16* c = reinterpret_cast<U16*>(&controller);
		stream.read(reinterpret_cast<char*>(c), sizeof(U16));
		return(stream);
	}

public:
	// 6 base values (4 integers + 2 floats)
	U16 signedness: 1,
		mixedStride: 1,
		type: 6,    // base typing (extra bits saved for future use)
		encoder: 5, // encoder bits (0 = uncompressed)
		uniform: 1, // triggered if all values are the same
		unused: 2;
};

}
}



#endif /* CORE_BASE_STREAMCONTAINERHEADERCONTROLLER_H_ */
