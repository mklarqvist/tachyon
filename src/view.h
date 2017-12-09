/*
Copyright (C) 2016-2017 Genome Research Ltd.
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
*/
#include <iostream>
#include <getopt.h>

#include "utility.h"
#include "core/TachyonReader.h"

void view_usage(void){
	programMessage();
	std::cerr <<
	"About:  Convert VCF/BCF->TWK/; subset and slice TWK/TWO data\n"
	"        Only biallelic diploid genotypes from SNVs will be retained\n"
	"Usage:  " << Tachyon::Constants::PROGRAM_NAME << " import [options] -i <in.vcf>/<in.bcf> -o <output.twk>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file prefix (required)\n"
	"  -c INT   checkpoint size in number of variants (default: 500)\n"
	"  -h FLOAT Hardy-Weinberg P-value cutoff (default: 0)\n"
	"  -m FLOAT Minor-genotype frequency (MGF) cutoff (default: 0)\n"
	"  -n FLOAT Missingness percentage cutoff (default: 0.2)\n"
	"  -s       Hide all program messages [null]\n";
}

int view(int argc, char** argv){
	if(argc == 1){
		programMessage();
		programHelpDetailed();
		return(1);
	}

	int c;
	if(argc < 2){
		import_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",		required_argument, 0,  'i' },
		{"output",		optional_argument, 0,  'o' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	SILENT = 0;

	while ((c = getopt_long(argc, argv, "i:o:c:s?", long_options, &option_index)) != -1){
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
		case 's':
			SILENT = 1;
			break;

		default:
			std::cerr << Tachyon::Helpers::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(input.length() == 0){
		std::cerr << Tachyon::Helpers::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << Tachyon::Helpers::timestamp("LOG") << "Calling view..." << std::endl;
	}

	Tachyon::Core::TachyonReader reader(input);
	//reader.settings.loadAll();

	reader.settings.loadMetaHot = true;
	//reader.settings.loadMetaCold = true;
	// Todo: deduplicate and move to function in settings class
	reader.settings.load_info_ID.push_back(5);
	reader.settings.load_info_ID.push_back(15);
	reader.settings.load_info_ID.push_back(16);
	reader.settings.load_info_ID.push_back(17);
	reader.settings.load_info_ID.push_back(18);
	reader.settings.load_info_ID.push_back(11);
	reader.settings.load_info_ID.push_back(24);
	reader.settings.load_info_ID.push_back(11);
	reader.settings.load_info_ID.push_back(25);

	reader.open(input);
	while(reader.nextBlock()){
		//reader.toVCF();
		reader.toVCFPartial();
		reader.block.clear();
	}

	return 0;
}
