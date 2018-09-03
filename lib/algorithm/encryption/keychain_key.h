#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_

#include "support/enums.h"

namespace tachyon {

struct KeychainKey {
public:
	KeychainKey() : field_id(0), encryption_type(YON_ENCRYPTION_NONE){}
	virtual ~KeychainKey(){}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const KeychainKey& key);
	friend std::ostream& operator<<(std::ostream& stream, const KeychainKey& key);
	friend std::istream& operator>>(std::istream& stream, KeychainKey& key);

	inline TACHYON_ENCRYPTION GetType(void) const{ return(TACHYON_ENCRYPTION(this->encryption_type)); }
	inline uint64_t& GetIdx(void){ return(this->field_id); }
	inline const uint64_t& GetIdx(void) const{ return(this->field_id); }

	virtual std::ostream& print(std::ostream& stream) =0;
	virtual KeychainKey* Clone() const =0;
	virtual io::BasicBuffer& AddToBuffer(io::BasicBuffer& buffer) const =0;
	virtual std::ostream& WriteToStream(std::ostream& stream) const =0;
	virtual std::istream& ReadFromStream(std::istream& stream) =0;

public:
	uint8_t   encryption_type;
	uint64_t  field_id;
};

template <uint8_t KeyLength = 32, uint8_t IVLength = 16>
struct KeychainKeyBase : public KeychainKey {
public:
	typedef KeychainKeyBase     self_type;
	typedef KeychainKey         parent_type;
	typedef io::BasicBuffer buffer_type;

public:
	KeychainKeyBase(){}
	virtual ~KeychainKeyBase(){}

	KeychainKeyBase(const KeychainKeyBase& other) :
		parent_type(other)
	{
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
	}

	KeychainKeyBase& operator=(const KeychainKeyBase& other){
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		return(*this);
	}

	virtual inline KeychainKey* Clone() const{ return new self_type(*this); }

	virtual io::BasicBuffer& AddToBuffer(io::BasicBuffer& buffer) const{
		buffer += this->encryption_type;
		buffer += this->field_id;
		buffer.Add((const char*)&this->key[0], KeyLength);
		buffer.Add((const char*)&this->iv[0],  IVLength);
		return(buffer);
	}

	virtual std::ostream& WriteToStream(std::ostream& stream) const{
		stream.write(reinterpret_cast<const char*>(&this->encryption_type), sizeof(uint8_t));
		stream.write(reinterpret_cast<const char*>(&this->field_id), sizeof(uint64_t));
		stream.write((const char*)&this->key[0], KeyLength);
		stream.write((const char*)&this->iv[0],  IVLength);
		return(stream);
	}

	virtual std::istream& ReadFromStream(std::istream& stream){
		//stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(uint8_t));
		stream.read(reinterpret_cast<char*>(&this->field_id), sizeof(uint64_t));
		stream.read((char*)&this->key[0], KeyLength);
		stream.read((char*)&this->iv[0],  IVLength);
		return(stream);
	}

	virtual std::ostream& print(std::ostream& stream){
		for(uint32_t i = 0; i < KeyLength; ++i) stream << std::hex << (int)this->key[i];
		stream << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) stream << std::hex << (int)this->iv[i];
		return(stream);
	}

	inline bool InitiateRandom(void){
		RAND_bytes(&this->key[0], KeyLength);
		RAND_bytes(&this->iv[0],  IVLength);
		return true;
	}

public:
	uint8_t key[KeyLength];  // 256 bit key
	uint8_t iv[IVLength];   // 128 bit initiation vector
};

template <uint8_t KeyLength = 32, uint8_t IVLength = 16, uint8_t TagLength = 16>
struct KeychainKeyGCM : public KeychainKeyBase<KeyLength, IVLength>{
public:
	typedef KeychainKeyGCM                   self_type;
	typedef KeychainKeyBase<KeyLength, IVLength> parent_type;

public:
	KeychainKeyGCM(){
		this->encryption_type = YON_ENCRYPTION_AES_256_GCM;
	}
	~KeychainKeyGCM(){}

	KeychainKeyGCM(const KeychainKeyGCM& other) :
		parent_type(other)
	{
		this->encryption_type = YON_ENCRYPTION_AES_256_GCM;
		memcpy(&this->tag[0], &other.tag[0], TagLength);
	}

	KeychainKeyGCM& operator=(const KeychainKeyGCM& other){
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		memcpy(&this->tag[0], &other.tag[0], TagLength);
		return(*this);
	}

	inline KeychainKey* Clone() const{ return new self_type(*this); }

	std::ostream& print(std::ostream& stream){
		for(uint32_t i = 0; i < KeyLength; ++i) stream << std::hex << (int)this->key[i];
		stream << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) stream << std::hex << (int)this->iv[i];
		stream << '\t';
		for(uint32_t i = 0; i < TagLength; ++i) stream << std::hex << (int)this->tag[i];
		stream << std::dec;
		return(stream);
	}

	io::BasicBuffer& AddToBuffer(io::BasicBuffer& buffer) const{
		buffer += this->encryption_type;
		buffer += this->field_id;
		buffer.Add((const char*)&this->key[0], KeyLength);
		buffer.Add((const char*)&this->iv[0],  IVLength);
		buffer.Add((const char*)&this->tag[0], TagLength);
		return(buffer);
	}

	std::ostream& WriteToStream(std::ostream& stream) const{
		stream.write(reinterpret_cast<const char*>(&this->encryption_type), sizeof(uint8_t));
		stream.write(reinterpret_cast<const char*>(&this->field_id), sizeof(uint64_t));
		stream.write((const char*)&this->key[0], KeyLength);
		stream.write((const char*)&this->iv[0],  IVLength);
		stream.write((const char*)&this->tag[0], TagLength);
		return(stream);
	}

	std::istream& ReadFromStream(std::istream& stream){
		//stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(uint8_t));
		stream.read(reinterpret_cast<char*>(&this->field_id), sizeof(uint64_t));
		stream.read((char*)&this->key[0], KeyLength);
		stream.read((char*)&this->iv[0],  IVLength);
		stream.read((char*)&this->tag[0], TagLength);
		return(stream);
	}

public:
	uint8_t tag[TagLength];  // validation tag for gcm
};

}



#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_ */
