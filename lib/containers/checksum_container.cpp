#include "checksum_container.h"

namespace tachyon{
namespace containers{

ChecksumContainer::ChecksumContainer(void) :
	n_entries(0),
	n_capacity(1000),
	entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::ChecksumContainer(const size_type capacity) :
	n_entries(0),
	n_capacity(capacity),
	entries_(static_cast<pointer>(::operator new[](this->capacity()*sizeof(value_type))))
{

}

ChecksumContainer::~ChecksumContainer() {
	for (std::size_t i = 0; i < this->n_entries; ++i)
		((this->entries_ + i)->~DigitalDigestPair)();

	::operator delete[](static_cast<void*>(this->entries_));
}

bool ChecksumContainer::Allocate(const size_type n_entries) {
	if (n_entries > this->capacity())
		return false;

	this->n_entries = n_entries;
	for (size_type i = 0; i < this->size(); ++i) {
		new( &this->entries_[this->n_entries] ) value_type( );
		this->entries_[i].compressed.initialize();
		this->entries_[i].uncompressed.initialize();
	}

	return true;
}

void ChecksumContainer::operator+=(const_reference value) {
	if (this->size() + 1 == this->capacity()) {
		// is full
		std::cerr << "is full: need resize" << std::endl;
		return;
	}

	// Invoke copy-ctor
	new( &this->entries_[this->n_entries] ) value_type( value );
	this->entries_[this->n_entries].compressed.initialize();
	this->entries_[this->n_entries].uncompressed.initialize();
	++this->n_entries;
}

void ChecksumContainer::Finalize(void) {
	if (this->size() == 0) return;
	for (size_type i = 0; i < this->size(); ++i)
		this->at(i).finalize();
}

bool ChecksumContainer::Update(const block_type& block, const header_type& header) {
	for (uint32_t i = 0; i < block.footer.n_info_streams; ++i) {
		//assert(mapTable[block.index_entry.info_offsets[i].key] < this->size());
		if (!(*this)[block.footer.info_offsets[i].data_header.global_key].uncompressed.update(block.info_containers[i].data_uncompressed, block.info_containers[i].strides_uncompressed, block.info_containers[i].header.data_header.HasMixedStride())) {
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}

		assert(block.footer.info_offsets[i].data_header.global_key < this->size());
		if (!(*this)[block.footer.info_offsets[i].data_header.global_key].compressed.update(block.info_containers[i].data, block.info_containers[i].strides, block.info_containers[i].header.data_header.HasMixedStride())) {
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}
	}

	for (uint32_t i = 0; i < block.footer.n_format_streams; ++i) {
		if (!(*this)[block.footer.format_offsets[i].data_header.global_key].uncompressed.update(block.format_containers[i].data_uncompressed, block.format_containers[i].strides_uncompressed, block.format_containers[i].header.data_header.HasMixedStride())) {
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}

		if (!(*this)[block.footer.format_offsets[i].data_header.global_key].compressed.update(block.format_containers[i].data, block.format_containers[i].strides, block.format_containers[i].header.data_header.HasMixedStride())) {
			std::cerr << utility::timestamp("ERROR","DIGEST") << "Failed to update digest..." << std::endl;
			return false;
		}
	}
	return true;
}

std::ostream& operator<<(std::ostream& out, const ChecksumContainer& container) {
	out.write((const char* const)reinterpret_cast<const size_t* const>(&container.n_entries), sizeof(size_t));
	for (int32_t i = 0; i < container.size(); ++i)
		out << container[i];

	return(out);
}

std::ifstream& operator>>(std::ifstream& stream, ChecksumContainer& container) {
	stream.read((char*)reinterpret_cast<size_t*>(&container.n_entries), sizeof(size_t));
	container.Allocate(container.n_entries);
	for (int i = 0; i < container.size(); ++i)
		stream >> container[i];

	return(stream);
}

}

}
