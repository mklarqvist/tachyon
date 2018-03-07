#ifndef READ_IMPORTER_H_
#define READ_IMPORTER_H_

#include "string"
#include "istream"
#include "iostream"

#include "../algorithm/OpenHashTable.h"
#include "../algorithm/compression/compression_container.h"
#include "../algorithm/timer.h"
#include "../algorithm/compression/libzpaq.h"
// Handle errors in libzpaq and elsewhere
void libzpaq::error(const char* msg) {
  if (strstr(msg, "ut of memory")) throw std::bad_alloc();
  throw std::runtime_error(msg);
}
#include "../containers/datacontainer.h"

namespace tachyon {

class In: public libzpaq::Reader {
public:
	In(io::BasicBuffer& buffer) : iterator_pos(0), buffer(buffer){}
	~In(){ }
	int get() {
		if(this->iterator_pos + 1 == this->buffer.size()) return(-1); // eof
		assert(this->buffer[this->iterator_pos] != -1);
		return(this->buffer[this->iterator_pos++]);
	}  // returns byte 0..255 or -1 at EOF

	void reset(void){
		this->buffer.reset();
		this->iterator_pos = 0;
	}

	size_t           iterator_pos;
	io::BasicBuffer& buffer;
 };

 class Out: public libzpaq::Writer {
 public:
	 Out(io::BasicBuffer& buffer) : buffer(buffer){}
	 ~Out(){ }
   void put(int c) { this->buffer += (BYTE)c; }  // writes 1 byte 0..255
   void reset(void){ this->buffer.reset(); }

   io::BasicBuffer& buffer;
 };

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
				container += (SBYTE)this->data[this->positions[j]][i];
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

		for(U32 i = 0; i < 6; ++i){
			if(i == 1){
				this->deltaEncode(this->nameContainer[i]);
			}
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
		container.header.uLength = container.buffer_data_uncompressed.size();

		// If we have mixed striding
		if(container.header.hasMixedStride()){
			// Reformat stream to use as small word size as possible
			container.reformatStride();
			container.header_stride.uLength = container.buffer_strides_uncompressed.size();
		}
	}

	void reset(void){
		qualContainer.reset();
		qualContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
		qualContainer.header.controller.signedness = true;

		basesContainer.reset();
		basesContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
		basesContainer.header.controller.signedness = true;

		for(U32 i = 0; i < 8; ++i){
			//std::cerr << nameContainer[i].getSizeUncompressed() << std::endl;
			nameContainer[i].reset();
			nameContainer[i].setType(TACHYON_CORE_TYPE::YON_TYPE_32B);
			nameContainer[i].header.controller.signedness = true;
		}
	}

	void deltaEncode(containers::DataContainer& container){
		if(container.size() == 0)
			return;

		// Recode integer types only
		if(!(container.header.controller.type == YON_TYPE_32B && container.header.controller.signedness == 1)){
			return;
		}

		if(container.header.controller.uniform == true)
			return;

		// At this point all integers are S32
		const S32* const dat  = reinterpret_cast<const S32* const>(container.buffer_data_uncompressed.buffer);
		container.buffer_data += dat[0];
		for(U32 j = 1; j < container.n_additions; ++j){
			container.buffer_data += dat[j] - dat[j-1];
			//std::cerr << "diff: " << dat[j] - dat[j-1] << std::endl;
		}
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

		//hash::HashTable<U64, U64> mers3(5024);
		//hash::HashTable<U64, U64> mers4(5024);

		std::string line;
		U64 count = 0;

		ReadBlock block;
		algorithm::CompressionManager compression_manager;
		//ReadSorter sorter(10000, 100);

		std::ofstream tester("/home/mklarqvist/Downloads/testqual.yon", std::ios::binary | std::ios::out);


		U32 read = 0;       U32 read_local = 0;
		U64 cost_qual = 0;  U64 cost_qual_raw = 0;
		U64 cost_bases = 0; U64 cost_bases_raw = 0;
		U64 cost_names = 0; U64 cost_names_raw = 0;
		//io::BasicBuffer test_buffer(block.basesContainer.buffer_data_uncompressed);
		//io::BasicBuffer out_buffer(block.basesContainer.buffer_data);
		//block.basesContainer.buffer_data_uncompressed.resize(10000*100*2);
		//block.basesContainer.buffer_data.resize(10000*100*2);
		In in(block.basesContainer.buffer_data_uncompressed);
		Out out(block.basesContainer.buffer_data);

		//block.qualContainer.buffer_data_uncompressed.resize(10000*100*2);
		//block.qualContainer.buffer_data.resize(10000*100*2);
		In in_qual(block.qualContainer.buffer_data_uncompressed);
		Out out_qual(block.qualContainer.buffer_data);

		algorithm::Timer timer; timer.Start();

		while(getline(std::cin, line)){
			if(count % 4 == 1){
				//std::cerr << line << std::endl;
				if(block.basesContainer.n_additions == 0){
					block.basesContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
					block.basesContainer.setStrideSize(line.size());
				}

				for(U32 i = 0; i < line.size(); ++i){
					//sorter[read_local][i] = line[i];
					block.basesContainer += line[i];
					//test_buffer+= (BYTE)line[i];
				}
				++block.basesContainer;

				if(!block.basesContainer.checkStrideSize(line.size()))
					block.basesContainer.triggerMixedStride();
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
					block.nameContainer[0].setStrideSize(first[0].size());
				}
				block.nameContainer[0].buffer_data_uncompressed.Add(&first[0][0], first[0].size());
				if(!block.nameContainer[0].checkStrideSize(first[0].size()))
					block.nameContainer[0].triggerMixedStride();

				block.nameContainer[0].addStride(first[0].size());
				++block.nameContainer[0];

				//std::cerr << "int: " << first[1] << std::endl;
				block.nameContainer[1] += atoi(&first[1][0]);
				++block.nameContainer[1];


				// split second into parts
				std::vector<std::string> end = utility::split(parts[1], '/');
				//std::cerr << end[0] << std::endl;
				//std::cerr << "int: " << end[1] << std::endl;
				block.nameContainer[2] += atoi(&end[1][0]);
				++block.nameContainer[2];

				// split end part
				std::vector<std::string> final = utility::split(end[0], ':');
				//std::cerr << "string: " << final[0] << std::endl;
				if(block.nameContainer[3].n_additions == 0){
					block.nameContainer[3].setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
					block.nameContainer[3].setStrideSize(final[0].size());
				}
				block.nameContainer[3].buffer_data_uncompressed.Add(&final[0][0], final[0].size());
				if(!block.nameContainer[3].checkStrideSize(final[0].size()))
					block.nameContainer[3].triggerMixedStride();

				block.nameContainer[3].addStride(final[0].size());
				++block.nameContainer[3];

				U32 target_container = 4;
				for(U32 i = 1; i < final.size(); ++i){
					block.nameContainer[target_container] += atoi(&final[i][0]);
					++block.nameContainer[target_container++];
				}
			}

			if(count % 4 == 3){
				if(block.qualContainer.n_additions == 0){
					block.qualContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_8B);
					block.qualContainer.setStrideSize(line.size());
				}

				for(U32 i = 0; i < line.size(); ++i){
					block.qualContainer += (BYTE)(line[i] - 33);
					++block.qualContainer;
					//qual_buffer += (BYTE)(line[i] - 33);
				}
				++read; ++read_local;

				if(read % 10000 == 0 && read != 0){
					block.updateContainers();
					libzpaq::compress(&in, &out, "x0.3ci1");
					//std::cerr << "bases: " << in.buffer.size() << "->" << out.buffer.size() << std::endl;
					libzpaq::compress(&in_qual, &out_qual, "x0.3ci1");
					//std::cerr << "quality: " << in_qual.buffer.size() << "->" << out_qual.buffer.size() << std::endl;
					compression_manager.zstd_codec.setCompressionLevel(6);
					for(U32 i = 0; i < 8; ++i)
						compression_manager.zstd_codec.encode(block.nameContainer[i]);

					tester.write(out.buffer.data(), out.buffer.size());
					tester.write(out_qual.buffer.data(), out_qual.buffer.size());
					for(U32 i = 0; i < 8; ++i){
						tester << block.nameContainer[i];
						cost_names += block.nameContainer[i].getObjectSize();
						cost_names_raw += block.nameContainer[i].buffer_data_uncompressed.size();
					}
					cost_bases += block.basesContainer.getObjectSize();
					cost_bases_raw += block.basesContainer.buffer_data_uncompressed.size();
					cost_qual += block.qualContainer.getObjectSize();
					cost_qual_raw += block.qualContainer.buffer_data_uncompressed.size();
					read_local = 0;
					const U64 total_cost = cost_bases + cost_qual + cost_names;
					const U64 total_cost_raw = cost_bases_raw + cost_qual_raw + cost_names_raw;
					std::cerr << utility::timestamp("PROGRESS") << std::setw(14) << utility::ToPrettyString(read) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost_raw) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost) << ' ' << std::setw(14) << (double)total_cost_raw/total_cost << "-fold" << ' ' << timer.ElapsedString() << std::endl;
					block.reset();
					in.reset(); out.reset();
					in_qual.reset(); out_qual.reset();
				}
			}

			++count;
			//std::cerr << line << std::endl;
		}

		if(read_local){
			block.updateContainers();
			libzpaq::compress(&in, &out, "x0.3ci1");
			//std::cerr << "bases: " << in.buffer.size() << "->" << out.buffer.size() << std::endl;
			libzpaq::compress(&in_qual, &out_qual, "x0.3ci1");
			//std::cerr << "quality: " << in_qual.buffer.size() << "->" << out_qual.buffer.size() << std::endl;
			compression_manager.zstd_codec.setCompressionLevel(6);
			for(U32 i = 0; i < 8; ++i)
				compression_manager.zstd_codec.encode(block.nameContainer[i]);

			tester.write(out.buffer.data(), out.buffer.size());
			tester.write(out_qual.buffer.data(), out_qual.buffer.size());
			for(U32 i = 0; i < 8; ++i){
				tester << block.nameContainer[i];
				cost_names += block.nameContainer[i].getObjectSize();
				cost_names_raw += block.nameContainer[i].buffer_data_uncompressed.size();
			}
			cost_bases += block.basesContainer.getObjectSize();
			cost_bases_raw += block.basesContainer.buffer_data_uncompressed.size();
			cost_qual += block.qualContainer.getObjectSize();
			cost_qual_raw += block.qualContainer.buffer_data_uncompressed.size();
			read_local = 0;
			const U64 total_cost = cost_bases + cost_qual + cost_names;
			const U64 total_cost_raw = cost_bases_raw + cost_qual_raw + cost_names_raw;
			std::cerr << utility::timestamp("PROGRESS") << std::setw(14) << utility::ToPrettyString(read) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost_raw) << ' ' << std::setw(14) << utility::toPrettyDiskString(total_cost) << ' ' << std::setw(14) << (double)total_cost_raw/total_cost << "-fold" << ' ' << timer.ElapsedString() << std::endl;
			block.reset();
			in.reset(); out.reset();
			in_qual.reset(); out_qual.reset();
		}



		std::cerr << utility::timestamp("LOG", "PROGRESS") << cost_names_raw << "->" << cost_names << '\t' << cost_bases_raw << "->" << cost_bases << '\t' << cost_qual_raw << "->" << cost_qual << std::endl;
		tester.flush();
		tester.close();

		return true;
	}

private:


};

}



#endif /* READ_IMPORTER_H_ */
