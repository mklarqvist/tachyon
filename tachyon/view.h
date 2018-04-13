/*
Copyright (C) 2017-2018 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk819@cam.ac.uk>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
==============================================================================*/
#include <iostream>
#include <getopt.h>

#include <regex>

#include "utility.h"
#include "variant_reader.h"

void view_usage(void){

}

int view(int argc, char** argv){
	if(argc <= 2){
		programMessage();
		programHelpDetailed();
		return(1);
	}

	int c;
	if(argc < 2){
		view_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",		required_argument, 0,  'i' },
		{"output",		optional_argument, 0,  'o' },
		{"keychain",	optional_argument, 0,  'k' },
		{"filter",	optional_argument, 0,  'f' },
		{"noHeader",	no_argument, 0,  'H' },
		{"onlyHeader",	no_argument, 0,  'h' },
		{"dropFormat",	no_argument, 0,  'G' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	std::string keychain_file;
	std::vector<std::string> load_strings;
	SILENT = 0;
	bool dropFormat = false;
	bool headerOnly = false;
	bool showHeader = true;

	while ((c = getopt_long(argc, argv, "i:o:k:f:GshH?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			input = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'k':
			keychain_file = std::string(optarg);
			break;
		case 'f':
			load_strings.push_back(std::string(optarg));
			break;
		case 'G':
			dropFormat = true;
			break;
		case 's':
			SILENT = 1;
			break;
		case 'h':
			headerOnly = true;
			break;
		case 'H':
			showHeader = false;
			break;

		default:
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling view..." << std::endl;
	}

	tachyon::VariantReader reader;

	// temp
	if(keychain_file.size()){
		std::ifstream keychain_reader(keychain_file, std::ios::binary | std::ios::in);
		if(!keychain_reader.good()){
			std::cerr << "failed to open: " << keychain_file << std::endl;
			return 1;
		}
		char header[tachyon::constants::FILE_HEADER_LENGTH];
		keychain_reader.read(&header[0], tachyon::constants::FILE_HEADER_LENGTH);
		if(strncmp(&header[0], &tachyon::constants::FILE_HEADER[0], tachyon::constants::FILE_HEADER_LENGTH) != 0){
			std::cerr << "ILLEGAL KEY FILE" << std::endl;
			return 1;
		}

		S32 major = 0, minor = 0, release = 0;
		keychain_reader.read(reinterpret_cast<char*>(&major),   sizeof(S32));
		keychain_reader.read(reinterpret_cast<char*>(&minor),   sizeof(S32));
		keychain_reader.read(reinterpret_cast<char*>(&release), sizeof(S32));

		keychain_reader >> reader.keychain;
		if(!keychain_reader.good()){
			std::cerr << "failed to parse keychain" << std::endl;
			return 1;
		}
	}

	if(!reader.open(input)){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	if(headerOnly){
		std::cout << reader.header.literals << std::endl;
		reader.header.writeVCFHeaderString(std::cout, true);
		return(0);
	}

	if(showHeader){
		reader.header.writeVCFHeaderString(std::cout, !dropFormat);
	}

	reader.getSettings().loadAll(true);
	if(load_strings.size()){
		/*
		 * Reserved:
		 * 1) CHROM
		 * 2) POS
		 * 3) ....
		 */
		reader.getSettings().loadAllINFO(false);
		std::regex field_identifier_regex("^[A-Z_]{1,}$");
		std::vector<std::string> INFO_FIELDS;
		for(U32 i = 0; i < load_strings.size(); ++i){
			if(strncasecmp(load_strings[i].data(), "INFO=", 5) == 0){
				std::vector<std::string> ind = tachyon::utility::split(load_strings[i].substr(5,load_strings.size()-5), ',');
				for(U32 j = 0; j < ind.size(); ++j){
					std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
					if(std::regex_match(ind[j], field_identifier_regex)){
						std::cerr << ind[j] << " GOOD: " << std::regex_match(ind[j], field_identifier_regex) << std::endl;
						INFO_FIELDS.push_back(ind[j]);
						reader.getSettings().loadINFO(ind[j]);
					} else {
						std::cerr << "Failed: " << ind[j] << " in string " << load_strings[i] << std::endl;
					}
				}
			} else {
				std::cerr << "Unknown pattern: " << load_strings[i] << std::endl;
			}
		}
	}

	if(dropFormat){
		reader.getSettings().loadGenotypes(false);
		reader.getSettings().loadPPA_    = false;
		reader.getSettings().loadFORMAT_ = false;
	}

	U64 n_variants = 0;
	tachyon::algorithm::Timer timer;
	timer.Start();

	U32 n_blocks = 0;
	while(reader.nextBlock()){;
		n_variants += reader.outputVCFBuffer();
		++n_blocks;
	}

	std::cerr << "Blocks: " << n_blocks << std::endl;
	std::cerr << "Variants: " << tachyon::utility::ToPrettyString(n_variants) << " genotypes: " << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t' << timer.ElapsedString() << '\t' << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count())) << std::endl;

	return 0;
}
