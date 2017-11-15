#include "../core/TomahawkImporter.h"
#include "../core/TomahawkImportWriter.h"

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
	//this->encoder->DetermineBitWidth();
	this->permutator.setSamples(this->header_->samples);

	this->writer_.setHeader(reader.header);
	if(!this->writer_.Open(this->outputPrefix)){
		std::cerr << Helpers::timestamp("ERROR", "WRITER") << "Failed to open writer..." << std::endl;
		return false;
	}

	///
	/// TODO
	/// temp

	this->sort_order_helper.previous_position = 0;
	this->sort_order_helper.contigID = nullptr;
	this->sort_order_helper.prevcontigID = 0;
	this->writer_.totempole_entry.contigID = 0;
	this->writer_.totempole_entry.minPosition = 0;

	//std::cerr << "PPA_conventional\tPPA_best\tPPA_byte\tPPA_u16\tPPA_u32\tPPA_u64\trle_conventional\trle_best\trle_byte\trle_u16\trle_u32\trle_u64\tfd_rle_best_ppa_best\tmemory_savings_rle_ppa\tfc_rle_conventional_ppa_best" << std::endl;
	while(true){
		if(!reader.getVariants(this->checkpoint_size))
			break;

		S32 contigID = reader[0].body->CHROM;
		this->sort_order_helper.previous_position = reader[0].body->POS;
		this->sort_order_helper.contigID = &contigID;
		this->sort_order_helper.prevcontigID = contigID;
		this->writer_.totempole_entry.contigID = contigID;
		this->writer_.totempole_entry.minPosition = reader[0].body->POS;

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


		++this->header_->getContig(contigID); // update block count for this contigID
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
	this->writer_.flush(this->permutator);
	this->writer_.WriteFinal();

	if(this->writer_.getVariantsWritten() == 0){
		std::cerr << Helpers::timestamp("ERROR","IMPORT") << "Did not import any variants..." << std::endl;
		return false;
	}

	if(!SILENT)
		std::cerr << Helpers::timestamp("LOG", "WRITER") << "Wrote: " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.getVariantsWritten()))
														 << " variants to " << Helpers::NumberThousandsSeparator(std::to_string(this->writer_.blocksWritten()))
														 << " blocks..." << std::endl;

	return(true);
}

bool TomahawkImporter::parseBCFLine(bcf_entry_type& line){
	if(this->sort_order_helper.prevcontigID != line.body->CHROM){
		if(line.body->CHROM < this->sort_order_helper.prevcontigID){
			std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "Contigs are not sorted (" << (*this->header_)[this->sort_order_helper.prevcontigID].name << " > " << (*this->header_)[line.body->CHROM].name << ")..." << std::endl;
			exit(1);
		}

		if(!SILENT)
			std::cerr << Helpers::timestamp("LOG", "IMPORT") << "Switch detected: " << this->header_->getContig(this->sort_order_helper.prevcontigID).name << "->" << this->header_->getContig(line.body->CHROM).name << "..." << std::endl;

		this->sort_order_helper.previous_position = 0;

		// Get new contig value from header
		// and flush out data
		++this->header_->getContig(line.body->CHROM);
		//this->writer_.flush(this->permutator);
	}

	// Assert position is in range
	if(line.body->POS + 1 > this->header_->getContig(line.body->CHROM).length){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << (*this->header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > reported max size of contig (" << (*this->header_)[line.body->CHROM].length << ")..." << std::endl;
		return false;
	}

	// Assert file is ordered
	if(line.body->POS < this->sort_order_helper.previous_position){
		std::cerr << Helpers::timestamp("ERROR", "IMPORT") << "File is not sorted by coordinates (" << (*this->header_)[line.body->CHROM].name << ':' << line.body->POS+1 << " > " << (*this->header_)[line.body->CHROM].name << ':' << this->sort_order_helper.previous_position << ")..." << std::endl;
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
