#ifndef READ_IMPORTER_H_
#define READ_IMPORTER_H_

#include "string"
#include "istream"
#include "iostream"

#include "../algorithm/OpenHashTable.h"
#include "../algorithm/compression/compression_container.h"
#include "../algorithm/timer.h"
#include "../containers/datacontainer.h"
#include "../core/header/read_header.h"
#include "../core/footer/read_footer.h"

namespace tachyon {

class ReadSorter{
public:
	ReadSorter(const U32 n_rows, const U32 l_reads) :
		n_rows(n_rows),
		l_reads(l_reads),
		data(new SBYTE*[n_rows]),
		positions(new U32[n_rows]),
		bins(new U32*[5])
	{
		for(U32 i = 0; i < n_rows; ++i)
			this->data[i] = new SBYTE[l_reads];

		for(U32 i = 0; i < 5; ++i)
			this->bins[i] = new U32[n_rows];

		for(U32 i = 0; i < n_rows; ++i)
			this->positions[i] = i;

		memset(&this->p_i[0], 0, 5*sizeof(U32));
	}

	~ReadSorter(){
		for(U32 i = 0; i < n_rows; ++i){
			delete [] this->data[i];
		}
		delete [] this->data;

		for(U32 i = 0; i < 5; ++i){
			delete [] this->bins[i];
		}
		delete [] this->bins;
		delete [] this->positions;
	}

	void clear(void){
		for(U32 i = 0; i < n_rows; ++i)
			this->positions[i] = i;

		memset(&this->p_i[0], 0, 5*sizeof(U32));
	}

	inline SBYTE* operator[](const U32 p){ return(this->data[p]); }

	void build(containers::DataContainer& container){
		for(S32 i = this->l_reads; i >= 0; --i){
			for(U32 j = 0; j < this->n_rows; ++j){
				//std::cerr << i << "/" << j << "\t" << this->positions[j] << '\t' << this->data[i][j] << std::endl;
				//this->bins[0][this->p_i[0]++] = this->positions[j];

				switch(this->data[this->positions[j]][i]){
				case('A'): this->bins[0][this->p_i[0]++] = this->positions[j]; break;
				case('T'): this->bins[1][this->p_i[1]++] = this->positions[j]; break;
				case('G'): this->bins[2][this->p_i[2]++] = this->positions[j]; break;
				case('C'): this->bins[3][this->p_i[3]++] = this->positions[j]; break;
				case('N'): this->bins[4][this->p_i[4]++] = this->positions[j]; break;
				default: std::cerr << "not possible: " << (char)this->data[j][i] << " at " << i << "," << j << std::endl; exit(1);
				}

			}

			U32 cum_pos = 0;
			for(U32 j = 0; j < 5; ++j){
				for(U32 p = 0; p < this->p_i[j]; ++p){
					this->positions[cum_pos + p] = this->bins[j][p];
				}
				cum_pos += this->p_i[j];
			}
			memset(&this->p_i[0], 0, 5*sizeof(U32));
			assert(cum_pos == this->n_rows);
		}

		std::cerr << "row: " << this->n_rows << "x" << this->l_reads << "->" << this->n_rows*this->l_reads << std::endl;
		container.buffer_data_uncompressed.reset();
		for(U32 j = 0; j < this->n_rows; ++j){
			for(U32 i = 0; i < this->l_reads; ++i){
				container.Add((SBYTE)this->data[this->positions[j]][i]);
			}
		}
		std::cerr << container.buffer_data_uncompressed.size() << std::endl;
	}

	U32     n_rows;
	U32     l_reads;
	SBYTE** data;
	U32*    positions;
	U32**   bins;
	U32     p_i[5];
};

class ReadBlock{
public:
	ReadBlock(){
		this->qualContainer.resize(200000);
		this->basesContainer.resize(200000);
		for(U32 i = 0; i < 8; ++i){
			this->nameContainer[i].resize(200000);
		}
	}
	~ReadBlock(){}

	void updateContainers(void){
		this->__updateContainer(this->basesContainer);
		this->__updateContainer(this->qualContainer);

		this->__updateContainer(this->nameContainer[0]);
		this->deltaEncode(this->nameContainer[1]);
		this->__updateContainer(this->nameContainer[1]);

		for(U32 i = 2; i < 6; ++i){
			this->__updateContainer(this->nameContainer[i]);
		}

		for(U32 i = 6; i < 8; ++i){
			this->deltaEncode(this->nameContainer[i]);
			this->__updateContainer(this->nameContainer[i]);
		}
	}

	void __updateContainer(containers::DataContainer& container){
		// Check if stream is uniform in content
		container.checkUniformity();
		// Reformat stream to use as small word size as possible
		container.reformat();

		// Set uncompressed length
		container.header.data_header.uLength = container.buffer_data_uncompressed.size();

		// If we have mixed striding
		if(container.header.data_header.hasMixedStride()){
			// Reformat stream to use as small word size as possible
			container.reformatStride();
			container.header.stride_header.uLength = container.buffer_strides_uncompressed.size();
		}
	}

	void reset(void){
		qualContainer.reset();
		qualContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
		qualContainer.header.data_header.controller.signedness = true;

		basesContainer.reset();
		basesContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
		basesContainer.header.data_header.controller.signedness = true;

		for(U32 i = 0; i < 8; ++i){
			//std::cerr << nameContainer[i].getSizeUncompressed() << std::endl;
			nameContainer[i].reset();
			nameContainer[i].setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
			nameContainer[i].header.data_header.controller.signedness = true;
		}
	}

	void deltaEncode(containers::DataContainer& container){
		if(container.size() == 0)
			return;

		// Recode integer types only
		if(!(container.header.data_header.controller.type == YON_TYPE_32B && container.header.data_header.controller.signedness == 1)){
			return;
		}

		if(container.header.data_header.controller.uniform == true)
			return;

		// At this point all integers are S32
		const S32* const dat  = reinterpret_cast<const S32* const>(container.buffer_data_uncompressed.buffer);

		// check for uniformity except first
		if(container.n_additions > 1){
			bool is_uniform_delta = true;
			const S32 first = dat[0];
			const S32 test_diff = dat[1] - dat[0];
			for(U32 i = 2; i < container.n_additions; ++i){
				if(dat[i] - dat[i - 1] != test_diff){
					is_uniform_delta = false;
					break;
				}
			}

			if(is_uniform_delta){
				container.n_entries   = 1;
				container.n_additions = 1;
				// Data pointers are updated in case there is no reformatting
				// see StreamContainer::reformat()
				container.buffer_data_uncompressed.n_chars = sizeof(S32);
				container.header.data_header.uLength                   = sizeof(S32);
				container.header.data_header.cLength                   = sizeof(S32);
				container.header.data_header.controller.uniform        = true;
				container.header.data_header.controller.mixedStride    = false;
				container.header.data_header.controller.encoder        = YON_ENCODE_NONE;


				return;
			}
		}

		container.buffer_data += dat[0];
		//std::cerr << dat[0];
		for(U32 j = 1; j < container.n_additions; ++j){
			container.buffer_data += dat[j] - dat[j-1];
			//std::cerr << ' ' << dat[j]-dat[j-1];
			//std::cerr << "diff: " << dat[j] - dat[j-1] << std::endl;
		}
		//std::cerr << std::endl;
		memcpy(container.buffer_data_uncompressed.data(),
				container.buffer_data.data(),
				container.buffer_data.size());

	}

public:
	containers::DataContainer basesContainer;
	containers::DataContainer qualContainer;
	containers::DataContainer nameContainer[8];
};

class ReadImporter{
private:
	typedef ReadImporter self_type;

public:
	ReadImporter(){}
	~ReadImporter(){}

	bool open(std::string file){
		//std::ifstream stream;
		//stream.open(file, std::ios::binary);
		//if(!stream.good()){
		//	std::cerr << "failed to open: " << file << std::endl;
		//}

		std::string line;
		U64 count = 0;

		ReadBlock block; block.reset();
		algorithm::CompressionManager compression_manager;
		//ReadSorter sorter(10000, 100);

		std::ofstream tester("/home/mklarqvist/Downloads/testqual.yon", std::ios::binary | std::ios::out);
		core::ReadHeader header;
		core::ReadFooter footer;
		header.write(tester);

		U32 read = 0;       U32 read_local = 0;
		U64 cost_qual = 0;  U64 cost_qual_raw = 0;
		U64 cost_bases = 0; U64 cost_bases_raw = 0;
		U64 cost_names = 0; U64 cost_names_raw = 0;

		//In in(block.basesContainer.buffer_data_uncompressed);
		//Out out(block.basesContainer.buffer_data);
		//In in_qual(block.qualContainer.buffer_data_uncompressed);
		//Out out_qual(block.qualContainer.buffer_data);

		algorithm::Timer timer; timer.Start();
		const U32 block_size = 50000;

		std::cerr << utility::timestamp("PROGRESS") << std::setw(14) << "Reads" << ' ' << std::setw(14) << "Input" << ' ' << std::setw(14) << "Output" << ' ' << std::setw(9+5) << "Compression" << ' ' << "Time" << std::endl;

		while(getline(std::cin, line)){
			if(count % 4 == 1){
				//std::cerr << line << std::endl;
				if(block.basesContainer.n_additions == 0){
					block.basesContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
				}

				for(U32 i = 0; i < line.size(); ++i){
					block.basesContainer.AddCharacter(line[i]);
				}
				++block.basesContainer;
				block.basesContainer.addStride(line.size());
			}

			if(count % 4 == 0){
				//std::cerr << line << std::endl;
				std::vector<std::string> parts = utility::split(line, ' ');
				//std::cerr << parts[0] << std::endl;
				//std::cerr << parts[1] << std::endl;
				// split first into parts
				std::vector<std::string> first = utility::split(parts[0], '.');

				//std::cerr << "string: " << first[0] << std::endl;
				if(block.nameContainer[0].n_additions == 0){
					block.nameContainer[0].setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
				}
				block.nameContainer[0].buffer_data_uncompressed.Add(&first[0][0], first[0].size());
				block.nameContainer[0].addStride(first[0].size());
				++block.nameContainer[0];

				//std::cerr << "int: " << first[1] << std::endl;
				block.nameContainer[1].Add(atoi(&first[1][0]));
				++block.nameContainer[1];


				// split second into parts
				std::vector<std::string> end = utility::split(parts[1], '/');
				//std::cerr << end[0] << std::endl;
				//std::cerr << "int: " << end[1] << std::endl;
				block.nameContainer[2].Add(atoi(&end[1][0]));
				++block.nameContainer[2];

				// split end part
				std::vector<std::string> final = utility::split(end[0], ':');
				//std::cerr << "string: " << final[0] << std::endl;
				if(block.nameContainer[3].n_additions == 0){
					block.nameContainer[3].setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
				}
				block.nameContainer[3].buffer_data_uncompressed.Add(&final[0][0], final[0].size());
				block.nameContainer[3].addStride(final[0].size());
				++block.nameContainer[3];

				U32 target_container = 4;
				for(U32 i = 1; i < final.size(); ++i){
					block.nameContainer[target_container].Add(atoi(&final[i][0]));
					++block.nameContainer[target_container++];
				}
			}

			if(count % 4 == 3){
				if(block.qualContainer.n_additions == 0){
					block.qualContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_8B);
				}

				for(U32 i = 0; i < line.size(); ++i){
					block.qualContainer.Add((BYTE)(line[i] - 33));
					++block.qualContainer;
				}
				++read; ++read_local;

				if(read % block_size == 0 && read != 0){
					for(U32 i = 0; i < 8; ++i)
						cost_names_raw += block.nameContainer[i].buffer_data_uncompressed.size();

					block.updateContainers();
					compression_manager.zpaq_codec.compress(block.basesContainer, false);
					compression_manager.zpaq_codec.compress(block.qualContainer, false);
					compression_manager.zstd_codec.setCompressionLevel(6);
					for(U32 i = 0; i < 8; ++i){
						if(i == 1) compression_manager.zpaq_codec.compress(block.nameContainer[i], false);
						else compression_manager.zstd_codec.compress(block.nameContainer[i]);
					}
					tester << block.basesContainer;
					tester << block.qualContainer;
					for(U32 i = 0; i < 8; ++i){
						tester << block.nameContainer[i];
						cost_names += block.nameContainer[i].getObjectSize();
					}
					cost_bases += block.basesContainer.getObjectSize();
					cost_bases_raw += block.basesContainer.buffer_data_uncompressed.size();
					cost_qual += block.qualContainer.getObjectSize();
					cost_qual_raw += block.qualContainer.buffer_data_uncompressed.size();
					const U64 total_cost = cost_bases + cost_qual + cost_names;
					const U64 total_cost_raw = cost_bases_raw + cost_qual_raw + cost_names_raw;
					std::cerr << utility::timestamp("PROGRESS") << std::setw(14) << utility::ToPrettyString(read) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost_raw) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost) << ' ' << std::setw(9) << (double)total_cost_raw/total_cost << "-fold" << ' ' << timer.ElapsedString() << std::endl;
					block.reset();

					++footer.n_blocks;
					footer.n_reads += block_size;
					tester.write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));
					read_local = 0;
				}
			}
			++count;
		}

		if(read_local){
			for(U32 i = 0; i < 8; ++i)
				cost_names_raw += block.nameContainer[i].buffer_data_uncompressed.size();

			block.updateContainers();
			compression_manager.zpaq_codec.compress(block.basesContainer, false);
			compression_manager.zpaq_codec.compress(block.qualContainer, false);
			compression_manager.zstd_codec.setCompressionLevel(6);
			for(U32 i = 0; i < 8; ++i)
				compression_manager.zstd_codec.compress(block.nameContainer[i]);

			tester << block.basesContainer;
			tester << block.qualContainer;
			for(U32 i = 0; i < 8; ++i){
				tester << block.nameContainer[i];
				cost_names += block.nameContainer[i].getObjectSize();
			}
			cost_bases += block.basesContainer.getObjectSize();
			cost_bases_raw += block.basesContainer.buffer_data_uncompressed.size();
			cost_qual += block.qualContainer.getObjectSize();
			cost_qual_raw += block.qualContainer.buffer_data_uncompressed.size();
			const U64 total_cost = cost_bases + cost_qual + cost_names;
			const U64 total_cost_raw = cost_bases_raw + cost_qual_raw + cost_names_raw;
			std::cerr << utility::timestamp("PROGRESS") << std::setw(14) << utility::ToPrettyString(read) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost_raw) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost) << ' ' << std::setw(9) << (double)total_cost_raw/total_cost << "-fold" << ' ' << timer.ElapsedString() << std::endl;
			block.reset();

			++footer.n_blocks;
			footer.n_reads += read_local;
			tester.write(reinterpret_cast<const char*>(&constants::TACHYON_BLOCK_EOF), sizeof(U64));
			read_local = 0;
		}
		footer.offset_end_of_data = tester.tellp();
		tester << footer;
		std::cerr << utility::timestamp("LOG", "FINAL") << utility::ToPrettyString(footer.n_blocks) << " blocks with " << utility::ToPrettyString(footer.n_reads) << " reads" << std::endl;
		std::cerr << utility::timestamp("LOG", "FINAL") << "Names:   " << cost_names_raw << "\t" << cost_names << '\t' << (float)cost_names_raw/cost_names << "-fold (" << (float)cost_names/cost_names_raw << ")"<< std::endl;
		std::cerr << utility::timestamp("LOG", "FINAL") << "Bases:   " << cost_bases_raw << "\t" << cost_bases << '\t' << (float)cost_bases_raw/cost_bases << "-fold (" << (float)cost_bases/cost_bases_raw << ")" << std::endl;
		std::cerr << utility::timestamp("LOG", "FINAL") << "Quality: " << cost_qual_raw << "\t" << cost_qual << '\t' << (float)cost_qual_raw/cost_qual << "-fold (" << (float)cost_qual/cost_qual_raw << ")" << std::endl;
		tester.flush();
		tester.close();
		return true;
	}

private:


};

}



#endif /* READ_IMPORTER_H_ */
