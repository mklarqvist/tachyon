#include "yon_tachyonheader.h"

namespace tachyon{
namespace core{

TachyonHeader::TachyonHeader(void) :
	version_major(0),
	version_minor(0),
	version_patch(0),
	n_contigs(0),
	n_samples(0),
	n_entries(0),
	contigs(nullptr),
	samples(nullptr),
	entries(nullptr),
	mapTable(nullptr),
	htable_contigs(nullptr),
	htable_samples(nullptr),
	htable_entries(nullptr)
{}

TachyonHeader::~TachyonHeader(){
	delete [] this->contigs;
	delete [] this->samples;
	delete [] this->entries;
	delete [] this->mapTable;
	delete this->htable_contigs;
	delete this->htable_samples;
	delete this->htable_entries;
}

const bool TachyonHeader::getContig(const std::string& p, contig_type*& target) const{
	if(this->htable_contigs == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_contigs->GetItem(&p[0], &p, ret, p.size())){
		target = &this->contigs[*ret];
		return true;
	}
	return false;
}

const bool TachyonHeader::getSample(const std::string& p, sample_type*& target) const{
	if(this->htable_samples == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_samples->GetItem(&p[0], &p, ret, p.size())){
		target = &this->samples[*ret];
		return true;
	}
	return false;
}

const bool TachyonHeader::getEntry(const std::string& p, map_entry_type*& target) const{
	if(this->htable_entries == nullptr) return false;
	S32* ret = nullptr;
	if(this->htable_entries->GetItem(&p[0], &p, ret, p.size())){
		target = &this->entries[*ret];
		return true;
	}
	return false;
}

bool TachyonHeader::buildMapTable(void){
	if(this->n_entries == 0)
		return false;

	delete [] this->mapTable;

	S32 largest_idx = -1;
	for(U32 i = 0; i < this->n_entries; ++i){
		if(this->entries[i].IDX > largest_idx)
			largest_idx = this->entries[i].IDX;
	}

	this->mapTable = new U32[largest_idx + 1];
	memset(this->mapTable, 0, sizeof(U32)*(largest_idx+1));
	S32 localID = 0;
	for(U32 i = 0; i < this->n_entries; ++i){
		this->mapTable[this->entries[i].IDX] = localID++;
		//std::cerr << i << "->" << this->mapTable[i] << std::endl;
	}

	return true;
}

bool TachyonHeader::buildHashTables(void){
	if(this->n_contigs){
		if(this->n_contigs*2 < 5012){
			this->htable_contigs = new hash_table_type(5012);
		} else
			this->htable_contigs = new hash_table_type(this->n_contigs*2);

		for(S32 i = 0; i < this->n_contigs; ++i){
			this->htable_contigs->SetItem(&this->contigs[i].name[0], &this->contigs[i].name, i, this->contigs[i].name.size());
		}
	}

	if(this->n_samples){
		if(this->n_samples*2 < 5012){
			this->htable_samples = new hash_table_type(5012);
		} else
			this->htable_samples = new hash_table_type(this->n_samples*2);

		for(S32 i = 0; i < this->n_samples; ++i){
			this->htable_samples->SetItem(&this->samples[i].name[0], &this->samples[i].name, i, this->samples[i].name.size());
		}
	}

	if(this->n_entries){
		if(this->n_entries*2 < 5012){
			this->htable_entries = new hash_table_type(5012);
		} else
			this->htable_entries = new hash_table_type(this->n_entries*2);

		for(S32 i = 0; i < this->n_entries; ++i){
			this->htable_entries->SetItem(&this->entries[i].ID[0], &this->entries[i].ID, i, this->entries[i].ID.size());
		}
	}

	return true;
}

}
}
