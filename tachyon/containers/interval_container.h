#ifndef CONTAINERS_INTERVAL_CONTAINER_H_
#define CONTAINERS_INTERVAL_CONTAINER_H_

#include <regex>

#include "../support/MagicConstants.h"
#include "../support/type_definitions.h"
#include "../support/helpers.h"
#include "../third_party/intervalTree.h"
#include "../index/index.h"
#include "../core/header/variant_header.h"
#include "../core/meta_entry.h"

namespace tachyon{
namespace containers{

class IntervalContainer {
private:
	typedef IntervalContainer      self_type;
    typedef std::size_t            size_type;
    typedef algorithm::Interval<U32, S64> interval_type;
    typedef algorithm::IntervalTree<U32, S64>  value_type;
    typedef value_type&            reference;
    typedef const value_type&      const_reference;
    typedef value_type*            pointer;
    typedef const value_type*      const_pointer;
    typedef index::Index           index_type;
    typedef index::IndexEntry      index_entry_type;
    typedef core::VariantHeader    header_type;
    typedef core::MetaEntry        meta_entry_type;

public:
    IntervalContainer() :
    	n_intervals_(0),
		n_entries_(0),
		__entries(nullptr)
	{

	}
    //IntervalContainer(const U32 n_contigs);
    //IntervalContainer(const U32 n_contigs, const std::vector<std::string>& interval_strings);

    ~IntervalContainer(void){
		for(std::size_t i = 0; i < this->n_entries_; ++i)
			(this->__entries + i)->~IntervalTree();

		::operator delete[](static_cast<void*>(this->__entries));
	}

	class iterator{
	private:
		typedef iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		reference operator*() const{ return *ptr_; }
		pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

	class const_iterator{
	private:
		typedef const_iterator self_type;
		typedef std::forward_iterator_tag iterator_category;

	public:
		const_iterator(pointer ptr) : ptr_(ptr) { }
		void operator++() { ptr_++; }
		void operator++(int junk) { ptr_++; }
		const_reference operator*() const{ return *ptr_; }
		const_pointer operator->() const{ return ptr_; }
		bool operator==(const self_type& rhs) const{ return ptr_ == rhs.ptr_; }
		bool operator!=(const self_type& rhs) const{ return ptr_ != rhs.ptr_; }
	private:
		pointer ptr_;
	};

    // Element access
    inline reference at(const size_type& position){ return(this->__entries[position]); }
    inline const_reference at(const size_type& position) const{ return(this->__entries[position]); }
    inline reference operator[](const size_type& position){ return(this->__entries[position]); }
    inline const_reference operator[](const size_type& position) const{ return(this->__entries[position]); }
    inline pointer data(void){ return(this->__entries); }
    inline const_pointer data(void) const{ return(this->__entries); }
    inline reference front(void){ return(this->__entries[0]); }
    inline const_reference front(void) const{ return(this->__entries[0]); }
    inline reference back(void){ return(this->__entries[this->n_entries_ - 1]); }
    inline const_reference back(void) const{ return(this->__entries[this->n_entries_ - 1]); }

    // Capacity
    inline const bool empty(void) const{ return(this->n_entries_ == 0); }
    inline const size_type& size(void) const{ return(this->n_entries_); }

    // Iterator
    inline iterator begin(){ return iterator(&this->__entries[0]); }
    inline iterator end(){ return iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator begin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator end() const{ return const_iterator(&this->__entries[this->n_entries_]); }
    inline const_iterator cbegin() const{ return const_iterator(&this->__entries[0]); }
    inline const_iterator cend() const{ return const_iterator(&this->__entries[this->n_entries_]); }

	inline const bool hasBlockList(void) const{ return(this->block_list_.size()); }
	inline std::vector<index_entry_type>& getBlockList(void){ return(this->block_list_); }
	inline const std::vector<index_entry_type>& getBlockList(void) const{ return(this->block_list_); }

    // Interpret
    bool validateIntervalStrings(std::vector<std::string>& interval_strings){
		if(interval_strings.size() == 0)
			return true;

		for(U32 i = 0; i < interval_strings.size(); ++i){
			// scrub whitespace
			//interval_strings[i].erase(remove_if(interval_strings[i].begin(), interval_strings[i].end(), isspace), interval_strings[i].end());
			interval_strings[i] = utility::remove_whitespace(interval_strings[i]);

			if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_ONLY )){
				//std::cerr << "chromosome onlu" << std::endl;
			} else if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_POSITION )){
				//std::cerr << "chromosome pos" << std::endl;
			} else if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_RANGE )){
				//std::cerr << "chromosome pos - pos" << std::endl;
			} else {
				std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << interval_strings[i] << std::endl;
				return false;
			}
		}

		return true;
	}

    /**<
	 * Parse interval strings. These strings have to match the regular expression
	 * patterns
	 * YON_REGEX_CONTIG_ONLY, YON_REGEX_CONTIG_POSITION, or YON_REGEX_CONTIG_RANGE
	 * @return Returns TRUE if successful or FALSE otherwise
	 */
	bool parseIntervals(std::vector<std::string>& interval_strings, const header_type& header, const index_type& index){
		// Intervals pass expression tests
		if(this->validateIntervalStrings(interval_strings) == false)
			return(false);

		// No intervals to parse
		if(interval_strings.size() == 0)
			return true;

		// Append given interval strings to internal vector of strings
		this->interval_strings_.insert( this->interval_strings_.end(), interval_strings.begin(), interval_strings.end() );

		// Assert that interval list data is of length n_contigs_;
		if(this->interval_list_.size() != header.getContigNumber())
			this->interval_list_ = std::vector< std::vector< interval_type > >(header.getContigNumber());


		// Parse each interval
		for(U32 i = 0; i < interval_strings.size(); ++i){
			interval_strings[i] = utility::remove_whitespace(interval_strings[i]);
			core::HeaderContig* contig = nullptr;

			// Chromosome only
			if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_ONLY )){
				std::cerr << "chromosome only" << std::endl;
				if(!header.getContig(interval_strings[i],contig)){
					std::cerr << "cant find contig: " << interval_strings[i] << std::endl;
					return(false);
				}

				std::cerr << "Parsed: " << interval_strings[i] << std::endl;
				this->interval_list_[contig->contigID].push_back(interval_type(0, contig->bp_length, contig->contigID));

			}
			// Chromosome:position
			else if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_POSITION )){
				std::cerr << "chromosome pos" << std::endl;
				std::vector<std::string> substrings = utility::split(interval_strings[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				U64 position = atof(substrings[1].data());
				std::cerr << "Parsed: " << substrings[0] << "," << position << std::endl;

				std::vector<index_entry_type> target_blocks = index.findOverlap(contig->contigID, position);
				this->block_list_.insert( this->block_list_.end(), target_blocks.begin(), target_blocks.end() );
				this->interval_list_[contig->contigID].push_back(interval_type(position, position, contig->contigID));


			}
			// Chromosome:position-position
			else if (std::regex_match (interval_strings[i], constants::YON_REGEX_CONTIG_RANGE )){
				std::cerr << "chromosome pos - pos" << std::endl;
				std::vector<std::string> substrings = utility::split(interval_strings[i], ':');
				if(substrings[0].size() == 0 || substrings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}

				if(!header.getContig(substrings[0],contig)){
					std::cerr << "cant find contig: " << substrings[0] << std::endl;
					return(false);
				}

				std::vector<std::string> position_strings = utility::split(substrings[1], '-');
				if(position_strings[0].size() == 0 || position_strings[1].size() == 0){
					std::cerr << "illegal form" << std::endl;
					return false;
				}
				U64 position_from = atof(position_strings[0].data());
				U64 position_to   = atof(position_strings[1].data());
				if(position_from > position_to) std::swap(position_from, position_to);

				std::cerr << "Parsed: " << substrings[0] << "," << position_from << "," << position_to << std::endl;

				std::vector<index_entry_type> target_blocks = index.findOverlap(contig->contigID, position_from, position_to);
				this->block_list_.insert( this->block_list_.end(), target_blocks.begin(), target_blocks.end() );
				this->interval_list_[contig->contigID].push_back(interval_type(position_from, position_to, contig->contigID));

			} else {
				std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << interval_strings[i] << std::endl;
				return false;
			}
			++this->n_intervals_;
		}

		if(this->block_list_.size() == 0)
			return true;

		return true;
	}

	bool build(const header_type& header){
		if(this->interval_list_.size() == 0) return true;

		// Dedupe blocks before building
		this->dedupeBlockList();

		this->n_entries_ = header.getContigNumber();
		this->__entries  = static_cast<pointer>(::operator new[](this->n_entries_*sizeof(value_type)));
		for(U32 i = 0; i < this->n_entries_; ++i){
			new( &this->__entries[i] ) value_type( std::move(this->interval_list_[i]) );
		}
		return true;
	}

	inline std::vector<interval_type> find_overlaps(const U32& contigID, const S64& start_position, const S64& end_position) const{
		if(contigID > this->size()) return(std::vector<interval_type>());
		return(this->at(contigID).findOverlapping(start_position, end_position));
	}

	inline std::vector<interval_type> find_overlaps(const meta_entry_type& meta_entry) const{
		if(meta_entry.contigID > this->size()) return(std::vector<interval_type>());
		return(this->at(meta_entry.contigID).findOverlapping(meta_entry.position, meta_entry.position + 1));
	}

private:
	void dedupeBlockList(void){
		if(this->block_list_.size() == 0) return;

		// Dedupe
		std::sort(this->block_list_.begin(), this->block_list_.end());
		std::vector<index_entry_type> block_list_deduped;
		block_list_deduped.push_back(this->block_list_[0]);
		for(U32 i = 1; i < this->block_list_.size(); ++i){
			if(block_list_deduped.back() != this->block_list_[i]){
				block_list_deduped.push_back(this->block_list_[i]);
			}
		}

		// Debug
		//std::cerr << "size: " << this->block_list__deduped.size() << std::endl;
		//for(U32 i = 0; i < this->block_list__deduped.size(); ++i){
		//	std::cerr << i << "\t";
		//	this->block_list__deduped[i].print(std::cerr);
		//	std::cerr << std::endl;
		//}

		// Assign deduped vector to settings reference
		this->block_list_ = block_list_deduped;
	}

private:
	U32 n_intervals_;
    std::vector<std::string>   interval_strings_;
    std::vector< std::vector<interval_type> > interval_list_;
    std::vector<index_entry_type> block_list_;
    size_t  n_entries_; // equal to number of contigs
    pointer __entries;  // interval trees
};

}
}

#endif /* CONTAINERS_INTERVAL_CONTAINER_H_ */
