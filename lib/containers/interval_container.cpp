#include "interval_container.h"

namespace tachyon{
namespace containers{

IntervalContainer::IntervalContainer() :
n_intervals_(0),
	n_entries_(0),
	__entries(nullptr)
{

}

IntervalContainer::~IntervalContainer(void){
	for(std::size_t i = 0; i < this->n_entries_; ++i)
		(this->__entries + i)->~IntervalTree();

	::operator delete[](static_cast<void*>(this->__entries));
}

// Interpret
bool IntervalContainer::validateIntervalStrings(std::vector<std::string>& interval_strings){
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

bool IntervalContainer::parseIntervals(std::vector<std::string>& interval_strings, const header_type& header, const index_type& index){
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

bool IntervalContainer::build(const header_type& header){
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

void IntervalContainer::dedupeBlockList(void){
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

}
}