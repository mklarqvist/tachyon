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
	programMessage(true);
	std::cerr <<
	"About:  Convert YON->VCF/BCF or custom output; provides subset and slice operators data\n"
	"Usage:  " << tachyon::constants::PROGRAM_NAME << " view [options] -i <in.yon>\n\n"
	"Options:\n"
	"  -i FILE   input YON file (required)\n"
	"  -o FILE   output file (- for stdout; default: -)\n"
	"  -k FILE   keychain with encryption keys (required if encrypted)\n"
	"  -f STRING interpreted filter string for slicing output (see manual)\n"
	"  -d CHAR   output delimiter (-c must be triggered)\n"
	"  -c        custom output format (ignores VCF/BCF specification rules)\n"
	"  -G        drop all FORMAT fields from output\n"
	"  -h/H      header only / no header\n"
	"  -s        Hide all program messages\n";
}

int view(int argc, char** argv){
	if(argc < 2){
		programMessage();
		programHelpDetailed();
		return(1);
	}

	int c;
	if(argc == 2){
		view_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",       required_argument, 0,  'i' },
		{"output",      optional_argument, 0,  'o' },
		{"keychain",    optional_argument, 0,  'k' },
		{"filter",      optional_argument, 0,  'f' },
		{"delimiter",   optional_argument, 0,  'd' },
		{"noHeader",    no_argument, 0,  'H' },
		{"onlyHeader",  no_argument, 0,  'h' },
		{"dropFormat",  no_argument, 0,  'G' },
		{"customFormat",no_argument, 0,  'c' },
		{"silent",      no_argument, 0,  's' },
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
	bool customDelimiter = false;
	char customDelimiterChar = 0;
	bool customOutputFormat = false;

	std::string temp;

	while ((c = getopt_long(argc, argv, "i:o:k:f:d:cGshH?", long_options, &option_index)) != -1){
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
		case 'c':
			customOutputFormat = true;
			break;
		case 'd':
			customDelimiter = true;
			temp = std::string(optarg);
			if(temp.size() != 1 && !(temp[0] == '\\' && temp.size() == 2)){
				std::cerr << "not a legal delimiter" << std::endl;
				return(1);
			}
			if(temp.size() == 1) customDelimiterChar = temp[0];
			else {
				if(temp[1] == 't') customDelimiterChar = '\t';
				else if(temp[1] == 'n') customDelimiterChar = '\n';
				else {
					std::cerr << "not a legal delimiter" << std::endl;
					return(1);
				}
			}
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

	if(load_strings.size()){
		reader.getSettings().custom_output_format = customOutputFormat;
		/*
		 * Reserved:
		 * 1) CHROM, CONTIG
		 * 2) POS, POSITION
		 * 3) REF, REFERENCE
		 * 4) ALT, ALTERNATIVE, ALTERNATE
		 * 5) QUAL, QUALITY
		 * 6) NAMES
		 */
		bool allGood = true;
		//reader.getSettings().loadAllINFO(false);
		//reader.getSettings().loadGenotypes(false);
		//reader.getSettings().load_ppa    = false;
		//reader.getSettings().load_format = false;

		std::regex field_identifier_regex("^[A-Z_0-9]{1,}$");
		for(U32 i = 0; i < load_strings.size(); ++i){
			std::vector<std::string> partitions = tachyon::utility::split(load_strings[i], ';');
			for(U32 p = 0; p < partitions.size(); ++p){
				partitions[p].erase(std::remove(partitions[p].begin(), partitions[p].end(), ' '), partitions[p].end()); // remove all spaces
				if(strncasecmp(partitions[p].data(), "INFO=", 5) == 0){
					std::vector<std::string> ind = tachyon::utility::split(partitions[p].substr(5,load_strings.size()-5), ',');
					for(U32 j = 0; j < ind.size(); ++j){
						//ind[j] = std::regex_replace(ind[j], std::regex("^ +| +$|( ) +"), "$1"); // remove excess white space
						std::transform(ind[j].begin(), ind[j].end(), ind[j].begin(), ::toupper); // transform to UPPERCASE
						if(std::regex_match(ind[j], field_identifier_regex)){
							const tachyon::core::HeaderMapEntry* map = reader.header.getInfoField(ind[j]);
							if(map == false){
								std::cerr << tachyon::utility::timestamp("ERROR") << "Failed: " << ind[j] << " in string " << partitions[p] << std::endl;
								allGood = false;
								continue;
							}
							reader.getSettings().loadINFO(ind[j]);
						} else {
							std::cerr << tachyon::utility::timestamp("ERROR") << "Failed: " << ind[j] << " in string " << partitions[p] << std::endl;
							allGood = false;
						}
					}
				} else if(strncasecmp(partitions[p].data(), "CONTIG", 6) == 0 && partitions[p].length() == 6){
					reader.getSettings().custom_output_controller.show_contig = true;
					reader.getSettings().load_contig = true;
				} else if(strncasecmp(partitions[p].data(), "POSITION", 8) == 0 && partitions[p].length() == 8){
					reader.getSettings().custom_output_controller.show_position = true;
					reader.getSettings().load_positons = true;
				} else if(strncasecmp(partitions[p].data(), "REF", 3) == 0 && partitions[p].length() == 3){
					reader.getSettings().custom_output_controller.show_ref = true;
					reader.getSettings().load_alleles = true;
				} else if(strncasecmp(partitions[p].data(), "ALT", 3) == 0 && partitions[p].length() == 3){
					reader.getSettings().custom_output_controller.show_alt = true;
					reader.getSettings().load_alleles = true;
				} else if(strncasecmp(partitions[p].data(), "QUALITY", 7) == 0 && partitions[p].length() == 7){
					reader.getSettings().custom_output_controller.show_quality = true;
					reader.getSettings().load_quality = true;
				} else if(strncasecmp(partitions[p].data(), "NAMES", 5) == 0 && partitions[p].length() == 5){
					reader.getSettings().custom_output_controller.show_names = true;
					reader.getSettings().load_names = true;
				} else {
					std::cerr << tachyon::utility::timestamp("ERROR") << "Unknown pattern: " << partitions[p] << std::endl;
					allGood = false;
				}
			}
		}

		if(allGood == false) return(1);
	} else {
		reader.getSettings().loadAll(true);

		if(dropFormat){
			reader.getSettings().loadGenotypes(false);
			reader.getSettings().load_ppa    = false;
			reader.getSettings().load_format = false;
		}
	}

	if(customDelimiter){
		if(customOutputFormat == false){
			std::cerr << "Have to trigger -c when using a custom separator" << std::endl;
			return(1);
		}
		reader.getSettings().setCustomDelimiter(customDelimiterChar);
	}

	tachyon::algorithm::Timer timer;
	timer.Start();

	if(showHeader){
		reader.getSettings().show_vcf_header = true;
	}

	U64 n_variants = 0;
	if(customOutputFormat) n_variants = reader.outputCustom();
	else n_variants = reader.outputVCF();

	//std::cerr << "Blocks: " << n_blocks << std::endl;
	std::cerr << "Variants: " << tachyon::utility::ToPrettyString(n_variants) << " genotypes: " << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t' << timer.ElapsedString() << '\t' << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count())) << std::endl;

	return 0;
}
