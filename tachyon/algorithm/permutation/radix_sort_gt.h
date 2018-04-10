#ifndef ALGORITHM_COMPRESSION_RADIXSORTGT_H_
#define ALGORITHM_COMPRESSION_RADIXSORTGT_H_

#include "../../io/bcf/BCFReader.h"
#include "../../algorithm/permutation/permutation_manager.h"
#include "../../core/genotype_summary.h"

namespace tachyon {
namespace algorithm {

/*
 * This class performs a radix sort on a
 * block of variant lines given they are
 * bi-allelic diploid.
 */
class RadixSortGT {
	typedef RadixSortGT        self_type;
	typedef bcf::BCFReader     bcf_reader_type;
	typedef bcf::BCFEntry      bcf_entry_type;
	typedef PermutationManager manager_type;

public:
	RadixSortGT();
	RadixSortGT(const U64 n_samples);
	~RadixSortGT();

	// Reset does NOT need to cast after each
	// iteration as values are overwritten
	// each cycle
	void reset(void);
	void setSamples(const U64 n_samples);

	// Construct given a reader with a block
	// of BCF entries loaded in it
	bool build(const bcf_reader_type& reader);
	bool update(const bcf_entry_type& entry);

	inline const U64& getSamples(void) const{ return(this->n_samples); }
	inline const U32& size(void) const{ return(this->position); }

	void test(const bcf_entry_type& entry, std::vector<containers::GenotypeSummary>& summaries){
		if(entry.hasGenotypes == false) return;

		if(entry.body->n_allele + 2 > summaries[0].n_alleles_){
			std::cerr << "too many alleles: " << entry.body->n_allele + 2 << "/" << (int)summaries[0].n_alleles_ << std::endl;
			return;
		}

		U32 internal_pos = entry.formatID[0].l_offset;
		U32 k = 0;
		for(U32 i = 0; i < 2*entry.body->n_sample; i += 2, ++k){
			const BYTE& fmt_type_value1 = *reinterpret_cast<const BYTE* const>(&entry.data[internal_pos++]);
			const BYTE& fmt_type_value2 = *reinterpret_cast<const BYTE* const>(&entry.data[internal_pos++]);
			BYTE alleleA = fmt_type_value2 >> 1;
			BYTE alleleB = fmt_type_value1 >> 1;
			alleleA += (alleleA != 0 ? 1 : 0);
			alleleB += (alleleB != 0 ? 1 : 0);

			++summaries[k].matrix_[alleleA][alleleB];
			//++this->vectorA_[alleleA];
			//++this->vectorB_[alleleB];
		}
	}

	void calculateDifference(const std::vector<containers::GenotypeSummary>& summaries){
		float** rows = new float*[this->n_samples];
		for(U32 i = 0; i < this->n_samples; ++i){
			rows[i] = new float[this->n_samples];
			memset(rows[i], 0, sizeof(float)*this->n_samples);
		}

		float mostSimilar = 0;
		U32 mostSimilar_i = 0, mostSimilar_j = 0;
		for(U32 i = 0; i < this->n_samples; ++i){ // individual i
			for(U32 j = i + 1; j < this->n_samples; ++j){ // individual j
				// compare summary statistics
				U64 similar = 0; U64 dissimilar = 0;
				for(U32 a = 0; a < 10; ++a){ // allele A
					for(U32 b = 0; b < 10; ++b){ // allele B
						if(summaries[i].matrix_[a][b] == summaries[j].matrix_[a][b]) ++similar;
						else ++dissimilar;
					}
				}
				rows[i][j] = (float)similar/(dissimilar+similar+1);
				if(rows[i][j] > mostSimilar){
					mostSimilar = rows[i][j];
					mostSimilar_i = i;
					mostSimilar_j = j;
				}
				//std::cerr << i << "/" << j << ":" << (float)similar/(dissimilar+similar+1) << std::endl;
			}
		}
		std::cerr << "Most similar: " << mostSimilar << "@" << mostSimilar_i <<"/" << mostSimilar_j << std::endl;

		bool* visited = new bool[this->n_samples];
		memset(visited, 0, sizeof(bool)*this->n_samples);
		U32* output = new U32[this->n_samples];

		// First point
		visited[mostSimilar_i] = true;
		visited[mostSimilar_j] = true;
		output[0] = mostSimilar_i;
		output[1] = mostSimilar_j;

		float nextSimilar = 0;
		U32 nextSimilarPerson = 0;
		for(U32 i = 2; i < this->n_samples; ++i){
			// From row
			for(U32 j = 0; j < i; ++j){
				if(rows[j][i] > nextSimilar && visited[j] == false){
					nextSimilar = rows[j][i];
					nextSimilarPerson = j;
				}
			}

			// From column
			for(U32 j = i+1; j < this->n_samples; ++j){
				if(rows[i][j] > nextSimilar && visited[j] == false){
					nextSimilar = rows[i][j];
					nextSimilarPerson = j;
				}
			}

			output[i] = nextSimilarPerson;
			visited[nextSimilarPerson] = true;
			nextSimilar = 0;
			nextSimilarPerson = 0;
		}

		for(U32 i = 0; i < this->n_samples; ++i){
			assert(visited[i]);
			(*this->manager)[i] = output[i];
			//std::cerr << (*this->manager)[i] << std::endl;
		}

		// temp dump
		/*
		for(U32 i = 0; i < this->n_samples; ++i){
			std::cout << rows[i][0];
			for(U32 j = 1; j < this->n_samples; ++j){
				std::cout << '\t' << rows[i][j];
			}
			std::cout << std::endl;
		}
		*/

		for(U32 i = 0; i < this->n_samples; ++i) delete rows[i];
		delete [] rows;
		delete [] visited;
		delete [] output;
		//exit(1);
	}

public:
	U64           n_samples; // total number of entries in file
	U32           position;  // number of entries parsed
	U32           p_i[9];    // number of entries in bin i
	BYTE*         GT_array;  // packed genotype array
	U32**         bins;      // bin i
	manager_type* manager;   // permutation manager
};

class GenotypeNearestNeighbour{
private:
	typedef GenotypeNearestNeighbour self_type;
	typedef bcf::BCFReader           bcf_reader_type;
	typedef bcf::BCFEntry            bcf_entry_type;

public:
	GenotypeNearestNeighbour(const U64& n_samples) :
		n_samples_(n_samples),
		output_permutation_order_(new U32[n_samples]),
		current_similarity_vector_(new U32[n_samples]),
		visited_samples_list_(new bool[n_samples]),
		n_matrix_rows_(0),
		n_biallelic_variants_(0),
		GT_matrix(new BYTE*[n_samples]),
		n_alt_count_(new U32[n_samples])
	{
		for(U32 i = 0; i < this->n_samples_; ++i)
			this->GT_matrix[i] = nullptr;

		this->reset();
	}

	~GenotypeNearestNeighbour(){
		delete [] this->output_permutation_order_;
		delete [] this->current_similarity_vector_;
		delete [] this->visited_samples_list_;
		delete [] this->n_alt_count_;
	}

	void resetIteration(void){
		memset(this->current_similarity_vector_, 0, sizeof(U32)*this->n_samples_);
	}

	void reset(void){
		memset(this->visited_samples_list_, 0, sizeof(bool)*this->n_samples_);
		memset(this->current_similarity_vector_, 0, sizeof(U32)*this->n_samples_);
		memset(this->output_permutation_order_, 0, sizeof(U32)*this->n_samples_);
		memset(this->n_alt_count_, 0, sizeof(U32)*this->n_samples_);
		this->n_biallelic_variants_ = 0;
	}

	void build(const bcf_reader_type& reader){
		if(reader.size() > this->n_matrix_rows_){
			for(U32 i = 0; i < this->n_samples_; ++i){
				delete [] this->GT_matrix[i];
				this->GT_matrix[i] = new BYTE[reader.size()];
			}
			this->n_matrix_rows_ = reader.size();
		}
		this->reset();

		// Construct matrix
		for(U32 i = 0; i < this->n_matrix_rows_; ++i){
			this->updateMatrix(reader[i]);
		}

		// Sort order
		std::vector< std::pair<U32,U32> > alt(this->n_samples_);
		for(U32 i = 0; i < this->n_samples_; ++i){
			alt[i].first = i;
			alt[i].second = this->n_alt_count_[i];
		}
		std::sort(alt.begin(), alt.end());

		// Greedy NN
		this->output_permutation_order_[0] = 0;
		this->visited_samples_list_[0] = true;
		for(U32 i = 1; i < this->n_samples_; ++i){
			this->update(reader[i], alt[i].first);
		}

		U32 cum_total = 0;
		for(U32 i = 0; i < this->n_samples_; ++i){
			cum_total += this->visited_samples_list_[i];
			std::cerr << this->output_permutation_order_[i] << ", ";
		}
		std::cerr << std::endl;
		assert(cum_total == this->n_samples_);
		// done
	}

	bool updateMatrix(const bcf_entry_type& entry){
		// Entry has to have genotypes
		if(entry.hasGenotypes == false)
			return false;

		// Has to be diploid and have no EOV (i.e. mixed ploidy)
		if(entry.gt_support.hasEOV || entry.gt_support.ploidy != 2)
			return false;

		// Has to be biallelic: otherwise skip
		if(!entry.isBiallelic())
			return false;

		if(entry.formatID[0].primitive_type != 1){
			std::cerr << "unexpected primitive: " << (int)entry.formatID[0].primitive_type << std::endl;
			return false;
		}

		U32 internal_pos = entry.formatID[0].l_offset;
		U32 k = 0;
		for(U32 i = 0; i < 2*this->n_samples_; i += 2, ++k){
			const SBYTE& fmt_type_value1 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
			const SBYTE& fmt_type_value2 = *reinterpret_cast<const SBYTE* const>(&entry.data[internal_pos++]);
			const BYTE packed = (bcf::BCF_UNPACK_GENOTYPE(fmt_type_value2) << 2) | bcf::BCF_UNPACK_GENOTYPE(fmt_type_value1);
			this->GT_matrix[k][this->n_biallelic_variants_] = packed;
			this->n_alt_count_ += (((fmt_type_value2 >> 1) - 1) == 0);
			this->n_alt_count_ += (((fmt_type_value1 >> 1) - 1) == 0);
		}
		++this->n_biallelic_variants_;
		return true;
	}

	void update(const bcf_entry_type& entry, const U32& reference_sample){
		this->resetIteration(); // reset local

		const BYTE* const reference = this->GT_matrix[reference_sample];
		for(U32 i = 0; i < this->n_samples_; ++i){
			if(i == reference_sample || this->visited_samples_list_[i]) continue;
			for(U32 r = 0; r < this->n_biallelic_variants_; ++r){
				if(reference[r] == this->GT_matrix[i][r])
					++current_similarity_vector_[i];
			}
		}

		S32 max = -1; S32 max_pos = -1;
		for(U32 i = 0; i < this->n_samples_; ++i){
			if((S32)current_similarity_vector_[i] > max && this->visited_samples_list_[i] == false){
				max = current_similarity_vector_[i];
				max_pos = i;

			}
			//std::cerr << current_similarity_vector_[i] << ", ";
		}
		//std::cerr << std::endl;
		//std::cerr << reference_sample << " -> " << "max: " << max << " at " << max_pos << std::endl;
		assert(max != -1);
		this->visited_samples_list_[(U32)max_pos] = true;
		this->output_permutation_order_[reference_sample] = max_pos;
	}

	const U32& operator[](const U32 position) const{ return(this->output_permutation_order_[position]); }
	const U32* const get(void) const{ return(this->output_permutation_order_); }

private:
	U64   n_samples_;
	U32*  output_permutation_order_;
	U32*  current_similarity_vector_;
	bool* visited_samples_list_;

	U32 n_matrix_rows_;
	U32 n_biallelic_variants_;
	BYTE** GT_matrix;
	U32*   n_alt_count_;
};

} /* namespace Algorithm */
} /* namespace Tomahawk */

#endif /* ALGORITHM_COMPRESSION_RADIXSORTGT_H_ */
