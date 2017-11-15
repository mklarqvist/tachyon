#include "TomahawkImporter.h"

namespace Tomahawk {

TomahawkImporter::TomahawkImporter(std::string inputFile, std::string outputPrefix, const U32 checkpoint) :
	checkpoint_size(checkpoint),
	block_flush_limit(65536),
	inputFile(inputFile),
	outputPrefix(outputPrefix),
	reader_(inputFile),
	writer_(),
	header_(nullptr),
	encoder(nullptr)
{}

TomahawkImporter::~TomahawkImporter(){
	delete this->encoder;
	// do not delete header: it might be a borrowed pointer
}

bool TomahawkImporter::Build(){
	std::ifstream temp(this->inputFile, std::ios::binary | std::ios::in);
	if(!temp.good()){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT")  << "Failed to open file..." << std::endl;
		return false;
	}
	char tempData[2];
	temp.read(&tempData[0], 2);
	temp.close();

	if((BYTE)tempData[0] == IO::Constants::GZIP_ID1 && (BYTE)tempData[1] == IO::Constants::GZIP_ID2){
		if(!this->BuildBCF()){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Failed build!" << std::endl;
			return false;
		}
	} else {
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Unknown file format!" << std::endl;
		return false;
	}
	return true;
}

bool TomahawkImporter::BuildBCF(void){
	bcf_reader_type reader;
	if(!reader.open(this->inputFile)){
		std::cerr << Helpers::timestamp("ERROR", "BCF")  << "Failed to open BCF file..." << std::endl;
		return false;
	}

	this->header_ = &reader.header;
	if(this->header_->samples == 0){
		std::cerr << Helpers::timestamp("ERROR", "BCF") << "No samples detected in header..." << std::endl;
		return false;
	}

	if(this->header_->samples == 1){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Cannot run " << Tomahawk::Constants::PROGRAM_NAME << " with a single sample..." << std::endl;
		return false;
	}

	// Spawn RLE controller
	this->encoder = new encoder_type(this->header_->samples);
	this->permutator.setSamples(this->header_->samples);

	this->writer_.setHeader(reader.header);
	if(!this->writer_.Open(this->outputPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	while(true){
		if(!reader.getVariants(this->checkpoint_size))
			break;

		// Reset permutate
		if(!this->permutator.build(reader)){
			std::cerr << "fail" << std::endl;
			return false;
		}

		for(U32 i = 0; i < reader.size(); ++i){
			if(!this->parseBCFLine(reader[i])){
				std::cerr << "failed to parse" << std::endl;
				return false;
			}
		}

		++this->header_->getContig(reader[0].body->CHROM); // update block count for this contigID
		this->writer_.flush(this->permutator);

		// Reset permutator
		this->permutator.reset();
	}

	// This only happens if there are no valid entries in the file
	if(this->sort_order_helper.contigID == nullptr){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	++this->header_->getContig(*this->sort_order_helper.contigID);
	//this->writer_.flush(this->permutator);
	//this->writer_.WriteFinal();

	/*
	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;
*/
	return(true);
}

bool TomahawkImporter::parseBCFLine(bcf_entry_type& line){
	// Assert position is in range
	if(line.body->POS + 1 > this->header_->getContig(line.body->CHROM).length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << (*this->header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > reported max size of contig (" << (*this->header_)[line.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Execute only if the line is simple (biallelic and SNP)
	if(line.isSimple()){

		// Flush if output block is over some size
		/*
		if(this->writer_.checkSize()){
			++this->header_->getContig(line.body->CHROM); // update block count for this contigID
			this->writer_.flush();

			this->writer_.TotempoleSwitch(line.body->CHROM, this->sort_order_helper.previous_position);
		}
		*/
		this->writer_.add(line, this->permutator.getPPA());
	}

	next:
	this->sort_order_helper.previous_position = line.body->POS;
	this->sort_order_helper.prevcontigID = line.body->CHROM;

	return true;
}

} /* namespace Tomahawk */
