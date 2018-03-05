#ifndef READ_IMPORTER_H_
#define READ_IMPORTER_H_

#include "string"
#include "istream"
#include "iostream"

#include "../algorithm/OpenHashTable.h"
#include "../algorithm/compression/compression_container.h"
#include "../containers/datacontainer.h"

namespace tachyon {

class ReadBlock{
public:
	ReadBlock(){
		this->qualContainer.resize(20000);
		this->basesContainer.resize(20000);
		for(U32 i = 0; i < 8; ++i){
			this->nameContainer[i].resize(200000);
		}
	}
	~ReadBlock(){}

	void updateContainers(void){
		this->__updateContainer(this->basesContainer);
		this->__updateContainer(this->qualContainer);
		for(U32 i = 0; i < 8; ++i)
			this->__updateContainer(this->nameContainer[i]);
	}

	void __updateContainer(containers::DataContainer& container){
		// Check if stream is uniform in content
		if(container.header.controller.type != tachyon::YON_TYPE_STRUCT){
			container.checkUniformity();
			// Reformat stream to use as small word size as possible
			container.reformat();
		}

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

		std::ofstream tester("/home/mklarqvist/Downloads/testqual.yon", std::ios::binary | std::ios::out);


		U32 read = 0;
		U64 cost_qual = 0;
		U64 cost_bases = 0;
		U64 cost_names = 0;
		while(getline(std::cin, line)){
			if(count % 4 == 1){
				//std::cerr << line << std::endl;
				if(block.basesContainer.n_additions == 0){
					block.basesContainer.setType(TACHYON_CORE_TYPE::YON_TYPE_CHAR);
					block.basesContainer.setStrideSize(line.size());
				}

				for(U32 i = 0; i < line.size(); ++i)
					block.basesContainer += line[i];

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
					//std::cerr << "int: " << final[i] << std::endl;
					block.nameContainer[target_container] += atoi(&final[i][0]);
					++block.nameContainer[target_container++];
				}
			}

			if(count % 4 == 3){
				for(U32 i = 0; i < line.size(); ++i){
					//quals[i] = line[i] - 33;
					//std::cout.write((char*)&quals[i], 1);
					block.qualContainer += (U32)(line[i] - 33);
					++block.qualContainer;
				}
				++read;

				if(read % 10000 == 0 && read != 0){
					block.updateContainers();
					compression_manager.zstd_codec.setCompressionLevel(19);
					compression_manager.zstd_codec.encode(block.qualContainer);
					tester << block.qualContainer;
					cost_qual += block.qualContainer.getObjectSize();
					std::cerr << utility::timestamp("LOG", "IMPORT") << read << ": " << block.qualContainer.getObjectSize() << " and " << block.qualContainer.buffer_data_uncompressed.size() << std::endl;

					compression_manager.zstd_codec.encode(block.basesContainer);
					tester << block.basesContainer;
					cost_bases += block.basesContainer.getObjectSize();
					std::cerr << utility::timestamp("LOG", "IMPORT") << "Bases: " << (int)block.basesContainer.getDataPrimitiveType() << ": " << block.basesContainer.getObjectSize() << " and " << block.basesContainer.buffer_data_uncompressed.size() << std::endl;

					for(U32 i = 0; i < 8; ++i){
						compression_manager.zstd_codec.encode(block.nameContainer[i]);
						tester << block.nameContainer[i];
						cost_names += block.nameContainer[i].getObjectSize();
						std::cerr << utility::timestamp("LOG", "IMPORT") << (int)block.nameContainer[i].getDataPrimitiveType() << ": " << block.nameContainer[i].getObjectSize() << " and " << block.nameContainer[i].buffer_data_uncompressed.size() << std::endl;
					}

					block.reset();
				}
			}

			++count;
			//std::cerr << line << std::endl;
		}

		if(block.qualContainer.getSizeUncompressed()){
			block.updateContainers();
			compression_manager.zstd_codec.setCompressionLevel(19);
			compression_manager.zstd_codec.encode(block.qualContainer);
			tester << block.qualContainer;
			cost_qual += block.qualContainer.getObjectSize();
			std::cerr << utility::timestamp("LOG", "IMPORT") << read << ": " << block.qualContainer.getObjectSize() << " and " << block.qualContainer.buffer_data_uncompressed.size() << std::endl;

			compression_manager.zstd_codec.encode(block.basesContainer);
			tester << block.basesContainer;
			cost_bases += block.basesContainer.getObjectSize();
			std::cerr << utility::timestamp("LOG", "IMPORT") << "Bases: " << (int)block.basesContainer.getDataPrimitiveType() << ": " << block.basesContainer.getObjectSize() << " and " << block.basesContainer.buffer_data_uncompressed.size() << std::endl;

			for(U32 i = 0; i < 8; ++i){
				compression_manager.zstd_codec.encode(block.nameContainer[i]);
				tester << block.nameContainer[i];
				cost_names += block.nameContainer[i].getObjectSize();
				std::cerr << utility::timestamp("LOG", "IMPORT") << (int)block.nameContainer[i].getDataPrimitiveType() << ": " << block.nameContainer[i].getObjectSize() << " and " << block.nameContainer[i].buffer_data_uncompressed.size() << std::endl;
			}

			block.reset();

		}


		return true;
	}

private:


};

}



#endif /* READ_IMPORTER_H_ */
