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
	"  -m        filtered data can match ANY number of requested fields\n"
	"  -M        filtered data must match ALL requested fields\n"
	"  -d CHAR   output delimiter (-c must be triggered)\n"
	"  -F STRING output format: can be either JSON,VCF,BCF, or CUSTOM (-c must be triggered)\n"
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
		{"filterAny",   no_argument, 0,  'm' },
		{"filterAll",   no_argument, 0,  'M' },
		{"delimiter",   optional_argument, 0,  'd' },
		{"output-type", optional_argument, 0,  'F' },
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
	bool filterAny = false;
	bool filterAll = false;

	std::string temp;

	while ((c = getopt_long(argc, argv, "i:o:k:f:d:F:cGshHmM?", long_options, &option_index)) != -1){
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

		case 'F':
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
			std::cerr << tachyon::utility::timestamp("ERROR") <<  "Failed to open keychain: " << keychain_file << "..." << std::endl;
			return 1;
		}

		keychain_reader >> reader.keychain;
		if(!keychain_reader.good()){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse keychain..." << std::endl;
			return 1;
		}
	}

	if(!reader.open(input)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << input << "..." << std::endl;
		return 1;
	}

	if(headerOnly){
		reader.header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		reader.header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
				  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + tachyon::utility::datetime();

		reader.header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE;

		std::cout << reader.header.literals << std::endl;
		reader.header.writeVCFHeaderString(std::cout, true);
		return(0);
	}

	// User provided '-f' string(s)
	if(load_strings.size()){
		if(!reader.getSettings().parseCommandString(load_strings, reader.header, customOutputFormat)){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse command..." << std::endl;
			return(1);
		}
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

	if(showHeader) reader.getSettings().show_vcf_header = true;
	else reader.getSettings().show_vcf_header = false;

	U64 n_variants = 0;
	if(customOutputFormat) n_variants = reader.outputCustom();
	else n_variants = reader.outputVCF();

	//std::cerr << "Blocks: " << n_blocks << std::endl;
	std::cerr << "Variants: " << tachyon::utility::ToPrettyString(n_variants) << " genotypes: " << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t' << timer.ElapsedString() << '\t' << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count())) << std::endl;

	return 0;
}
