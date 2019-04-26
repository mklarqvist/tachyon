#include "interval_container.h"

namespace tachyon{
namespace containers{

IntervalContainer::IntervalContainer() :
n_intervals_(0),
	n_entries_(0),
	__entries(nullptr)
{

}

IntervalContainer::~IntervalContainer(void) {
	for (std::size_t i = 0; i < this->n_entries_; ++i)
		(this->__entries + i)->~IntervalTree();

	::operator delete[](static_cast<void*>(this->__entries));
}

IntervalContainer::IntervalContainer(const self_type& other) :
	n_intervals_(other.n_intervals_),
	n_entries_(other.n_entries_),
    interval_strings_(other.interval_strings_),
    interval_list_(other.interval_list_),
    block_list_(other.block_list_),
	__entries(static_cast<pointer>(::operator new[](this->n_entries_*sizeof(value_type))))
{
	for (uint32_t i = 0; i < this->size(); ++i)
		new( &this->__entries[i] ) value_type( other.__entries[i] );
}

IntervalContainer::IntervalContainer(self_type&& other) noexcept :
	n_intervals_(other.n_intervals_),
	n_entries_(other.n_entries_),
	interval_strings_(std::move(other.interval_strings_)),
	interval_list_(std::move(other.interval_list_)),
	block_list_(std::move(other.block_list_)),
	__entries(nullptr)
{
	std::swap(this->__entries, other.__entries);
}

IntervalContainer& IntervalContainer::operator=(const self_type& other) {
	for (std::size_t i = 0; i < this->n_entries_; ++i)
		(this->__entries + i)->~IntervalTree();

	::operator delete[](static_cast<void*>(this->__entries));
	*this = IntervalContainer(other);
	return(*this);
}

IntervalContainer& IntervalContainer::operator=(self_type&& other) noexcept{
	for (std::size_t i = 0; i < this->n_entries_; ++i)
		(this->__entries + i)->~IntervalTree();

	::operator delete[](static_cast<void*>(this->__entries));
	this->__entries = nullptr;
	*this = IntervalContainer(std::move(other));
	return(*this);
}

// Interpret
bool IntervalContainer::ValidateIntervalStrings(std::vector<std::string>& interval_strings) {
	if (interval_strings.size() == 0)
		return true;

	for (uint32_t i = 0; i < interval_strings.size(); ++i) {
		// scrub whitespace
		//interval_strings[i].erase(remove_if (interval_strings[i].begin(), interval_strings[i].end(), isspace), interval_strings[i].end());
		interval_strings[i] = utility::remove_whitespace(interval_strings[i]);

		if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_ONLY )) {
			//std::cerr << "chromosome onlu" << std::endl;
		} else if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_POSITION )) {
			//std::cerr << "chromosome pos" << std::endl;
		} else if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_RANGE )) {
			//std::cerr << "chromosome pos - pos" << std::endl;
		} else {
			std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << interval_strings[i] << std::endl;
			return false;
		}
	}

	return true;
}

bool IntervalContainer::ParseIntervals(std::vector<std::string>& interval_strings, const header_type& header, const index_type& index) {
	// Intervals pass expression tests
	// No intervals to parse
	if (interval_strings.size() == 0)
		return(true);

	if (this->ValidateIntervalStrings(interval_strings) == false)
		return(false);

	// Append given interval strings to internal vector of strings
	this->interval_strings_.insert( this->interval_strings_.end(), interval_strings.begin(), interval_strings.end() );

	// Assert that interval list data is of length n_contigs_
	// Note that this will truncate previous entries if resizing occurs
	if (this->interval_list_.size() != header.GetNumberContigs())
		this->interval_list_ = std::vector< std::vector< interval_type > >(header.GetNumberContigs());

	// Parse each interval
	for (uint32_t i = 0; i < interval_strings.size(); ++i) {
		interval_strings[i] = utility::remove_whitespace(interval_strings[i]);
		const YonContig* contig = nullptr;

		// Chromosome only
		if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_ONLY )) {
			//std::cerr << "chromosome only" << std::endl;
			contig = header.GetContig(interval_strings[i]);
			if (contig == nullptr) {
				std::cerr << utility::timestamp("ERROR") << "Contig not defined in file: " << interval_strings[i] << std::endl;
				return(false);
			}

			//std::cerr << "Parsed: " << interval_strings[i] << " -> " << contig->name << ":" << contig->idx << std::endl;
			std::vector<index_entry_type> target_blocks = index.FindOverlap(contig->idx);
			this->block_list_.insert( this->block_list_.end(), target_blocks.begin(), target_blocks.end() );
			this->interval_list_[contig->idx].push_back(interval_type(0, contig->n_bases, contig->idx));
		}
		// Chromosome:position
		else if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_POSITION )) {
			//std::cerr << "chromosome pos" << std::endl;
			std::vector<std::string> substrings = utility::split(interval_strings[i], ':');
			if (substrings[0].size() == 0 || substrings[1].size() == 0) {
				std::cerr << utility::timestamp("ERROR") << "Illegal form: " << interval_strings[i] << std::endl;
				return false;
			}

			contig = header.GetContig(substrings[0]);
			if (contig == nullptr) {
				std::cerr << utility::timestamp("ERROR") << "Contig not defined in file: " << interval_strings[i] << std::endl;
				return(false);
			}

			uint64_t position = atof(substrings[1].data());
			//std::cerr << "Parsed: " << substrings[0] << "," << position << std::endl;

			std::vector<index_entry_type> target_blocks = index.FindOverlap(contig->idx, position);
			//std::cerr << "overlaps: " << target_blocks.size() << std::endl;
			this->block_list_.insert( this->block_list_.end(), target_blocks.begin(), target_blocks.end() );
			this->interval_list_[contig->idx].push_back(interval_type(position, position, contig->idx));
		}
		// Chromosome:position-position
		else if (std::regex_match (interval_strings[i], YON_REGEX_CONTIG_RANGE )) {
			//std::cerr << "chromosome pos - pos" << std::endl;
			std::vector<std::string> substrings = utility::split(interval_strings[i], ':');
			if (substrings[0].size() == 0 || substrings[1].size() == 0) {
				std::cerr << utility::timestamp("ERROR") << "Illegal form: " << std::endl;
				return false;
			}

			contig = header.GetContig(substrings[0]);
			if (contig == nullptr) {
				std::cerr << utility::timestamp("ERROR") << "Contig not defined in file: " << substrings[0] << std::endl;
				return(false);
			}

			std::vector<std::string> position_strings = utility::split(substrings[1], '-');
			if (position_strings[0].size() == 0 || position_strings[1].size() == 0) {
				std::cerr << utility::timestamp("ERROR") << "Illegal form: " << std::endl;
				return false;
			}
			uint64_t position_from = atof(position_strings[0].data());
			uint64_t position_to   = atof(position_strings[1].data());
			if (position_from > position_to) std::swap(position_from, position_to);

			//std::cerr << "Parsed: " << substrings[0] << "," << position_from << "," << position_to << std::endl;

			std::vector<index_entry_type> target_blocks = index.FindOverlap(contig->idx, position_from, position_to);
			this->block_list_.insert( this->block_list_.end(), target_blocks.begin(), target_blocks.end() );
			this->interval_list_[contig->idx].push_back(interval_type(position_from, position_to, contig->idx));

		} else {
			std::cerr << utility::timestamp("ERROR") << "Uninterpretable interval string: " << interval_strings[i] << std::endl;
			return false;
		}
		++this->n_intervals_;
	}

	if (this->block_list_.size() == 0)
		return(false);

	// Ascertain data is sorted according to YON blockID
	// This has the implication that the traversal of a sliced sorted file
	// will be sorted as output
	std::sort(this->block_list_.begin(), this->block_list_.end());

	return true;
}

bool IntervalContainer::Build(const header_type& header) {
	if (this->interval_list_.size() == 0) return true;

	// Dedupe blocks before building
	this->DedupeBlockList();

	this->n_entries_ = header.GetNumberContigs();
	this->__entries  = static_cast<pointer>(::operator new[](this->n_entries_*sizeof(value_type)));
	for (uint32_t i = 0; i < this->n_entries_; ++i) {
		new( &this->__entries[i] ) value_type( std::move(this->interval_list_[i]) );
	}
	return true;
}

void IntervalContainer::DedupeBlockList(void) {
	if (this->block_list_.size() == 0) return;

	// Dedupe
	std::sort(this->block_list_.begin(), this->block_list_.end());
	std::vector<index_entry_type> block_list_deduped;
	block_list_deduped.push_back(this->block_list_[0]);
	for (uint32_t i = 1; i < this->block_list_.size(); ++i) {
		if (block_list_deduped.back() != this->block_list_[i]) {
			block_list_deduped.push_back(this->block_list_[i]);
		}
	}

	// Debug
	//std::cerr << "size: " << this->block_list__deduped.size() << std::endl;
	//for (uint32_t i = 0; i < this->block_list__deduped.size(); ++i) {
	//	std::cerr << i << "\t";
	//	this->block_list__deduped[i].print(std::cerr);
	//	std::cerr << std::endl;
	//}

	// Assign deduped vector to settings reference
	this->block_list_ = block_list_deduped;
}

}
}
