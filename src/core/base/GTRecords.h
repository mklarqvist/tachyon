#ifndef CORE_BASE_GTRECORDS_H_
#define CORE_BASE_GTRECORDS_H_

namespace Tachyon{
namespace Core{

// We CANNOT place phasing template parameter
// since if we set it to 0 then we have a
// 0 width bit field which is illegal
// To solve this we introduce the TachyonRunNoPhase
// data structure below
template <class T, BYTE missing = 1>
struct __attribute__((packed)) TachyonRun{
private:
	typedef TachyonRun self_type;

public:
	TachyonRun();	// Disallowed ctor
	~TachyonRun(); // Disallowed dtor

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << (U32)entry.alleleA
			<< (entry.phasing ? '|' : '/')
			<< (U32)entry.alleleB;
		return(out);
	}

	T phasing: 1,
	  alleleA: (1 + missing),
	  alleleB: (1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing) + 1);
};

template <class T, BYTE missing = 1>
struct __attribute__((packed)) TachyonRunNoPhase{
private:
	typedef TachyonRunNoPhase self_type;

public:
	TachyonRunNoPhase();  // Disallowed ctor
	~TachyonRunNoPhase(); // Disallowed dtor

	T alleleA: (1 + missing),
	  alleleB: (1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing));
};

template <class T, BYTE missing = 1>
struct __attribute__((packed)) TachyonRunPacked{
private:
	typedef TachyonRunPacked self_type;

public:
	TachyonRunPacked();  // Disallowed ctor
	~TachyonRunPacked(); // Disallowed dtor

	T phasing: 1,
	  alleles: 2*(1 + missing),
	  runs:    sizeof(T)*8 - (2 * (1 + missing) + 1);
};

template <class T>
struct __attribute__((packed)) TachyonRunSimple{
private:
	typedef TachyonRunSimple<T> self_type;

public:
	TachyonRunSimple();  // Disallowed ctor
	~TachyonRunSimple(); // Disallowed dtor

	friend std::ostream& operator<<(std::ostream& out, const self_type& entry){
		out << (entry.alleleA == 0 ? '.' : (U16)entry.alleleA - 1)
			<< (entry.phasing ? '|' : '/')
			<< (entry.alleleB == 0 ? '.' : (U16)entry.alleleB - 1);
		return(out);
	}

	T phasing: 1,
	  alleleA: (sizeof(T)*8)/2 - 1,
	  alleleB: (sizeof(T)*8)/2 - 1,
	  unused:  1;
};

}
}

#endif /* CORE_BASE_GTRECORDS_H_ */
