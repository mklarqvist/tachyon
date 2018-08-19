#ifndef ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_
#define ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_

namespace tachyon{
namespace encryption{

template <uint8_t KeyLength = 32, uint8_t IVLength = 16>
struct KeychainKey{
public:
	typedef KeychainKey     self_type;
	typedef io::BasicBuffer buffer_type;

public:
	KeychainKey() : field_id(0), encryption_type(YON_ENCRYPTION_NONE){}
	virtual ~KeychainKey(){}

	KeychainKey(const KeychainKey& other) :
		field_id(other.field_id),
		encryption_type(other.encryption_type)
	{
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
	}

	KeychainKey& operator=(const KeychainKey& other){
		this->field_id = other.field_id;
		this->encryption_type = other.encryption_type;
		memcpy(&this->key[0], &other.key[0], KeyLength);
		memcpy(&this->iv[0],  &other.iv[0],  IVLength);
		return(*this);
	}

	virtual void print(void){
		for(uint32_t i = 0; i < KeyLength; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) std::cerr << std::hex << (int)this->iv[i];
	}

	inline bool initiateRandom(void){
		RAND_bytes(&this->key[0], KeyLength);
		RAND_bytes(&this->iv[0],  IVLength);
		return true;
	}

private:
	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.field_id;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], KeyLength);
		buffer.Add((const char*)&key.iv[0],  IVLength);
		return(buffer);
	}

public:
	uint64_t  field_id;
	uint8_t encryption_type;
	uint8_t key[KeyLength];  // 256 bit key
	uint8_t iv[IVLength];   // 128 bit initiation vector
};

template <uint8_t KeyLength = 32, uint8_t IVLength = 16, uint8_t TagLength = 16>
struct KeychainKeyGCM : public KeychainKey<KeyLength, IVLength>{
public:
	typedef KeychainKeyGCM                   self_type;
	typedef KeychainKey<KeyLength, IVLength> parent_type;

public:
	KeychainKeyGCM(){}
	~KeychainKeyGCM(){}

	KeychainKeyGCM(const KeychainKeyGCM& other) :
		parent_type(other)
	{
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

	void print(void){
		for(uint32_t i = 0; i < KeyLength; ++i) std::cerr << std::hex << (int)this->key[i];
		std::cerr << '\t';
		for(uint32_t i = 0; i < IVLength; ++i) std::cerr << std::hex << (int)this->iv[i];
		std::cerr << '\t';
		for(uint32_t i = 0; i < TagLength; ++i) std::cerr << std::hex << (int)this->tag[i];
		std::cerr << std::dec;
	}

private:
	friend std::ostream& operator<<(std::ostream& stream, const self_type& key){
		stream.write(reinterpret_cast<const char*>(&key.field_id), sizeof(uint64_t));
		stream.write(reinterpret_cast<const char*>(&key.encryption_type), sizeof(uint8_t));
		stream.write((const char*)&key.key[0], KeyLength);
		stream.write((const char*)&key.iv[0],  IVLength);
		stream.write((const char*)&key.tag[0], TagLength);
		return(stream);
	}

	friend std::istream& operator>>(std::istream& stream, self_type& key){
		stream.read(reinterpret_cast<char*>(&key.field_id), sizeof(uint64_t));
		stream.read(reinterpret_cast<char*>(&key.encryption_type), sizeof(uint8_t));
		stream.read((char*)&key.key[0], KeyLength);
		stream.read((char*)&key.iv[0],  IVLength);
		stream.read((char*)&key.tag[0], TagLength);
		return(stream);
	}

	friend io::BasicBuffer& operator+=(io::BasicBuffer& buffer, const self_type& key){
		buffer += key.field_id;
		buffer += key.encryption_type;
		buffer.Add((const char*)&key.key[0], KeyLength);
		buffer.Add((const char*)&key.iv[0],  IVLength);
		buffer.Add((const char*)&key.tag[0], TagLength);
		return(buffer);
	}

	friend io::BasicBuffer& operator>>(io::BasicBuffer& buffer, self_type& key){
		buffer >> key.field_id;
		buffer >> key.encryption_type;
		buffer.read((char*)&key.key[0], KeyLength);
		buffer.read((char*)&key.iv[0],  IVLength);
		buffer.read((char*)&key.tag[0], TagLength);
		return(buffer);
	}

public:
	uint8_t tag[TagLength];  // validation tag for gcm
};

}
}



#endif /* ALGORITHM_ENCRYPTION_KEYCHAIN_KEY_H_ */
