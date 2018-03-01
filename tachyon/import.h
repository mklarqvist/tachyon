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

#include "core/variant_importer.h"
#include "utility.h"

void import_usage(void){
	programMessage();
	std::cerr <<
	"Brief:  Convert VCF/BCF -> YON/\n"
	"Usage:  " << tachyon::constants::PROGRAM_NAME << " import [options] -i <in.bcf> -o <output.yon>\n\n"
	"Options:\n"
	"  -i FILE  input Tomahawk (required)\n"
	"  -o FILE  output file prefix (required)\n"
	"  -c INT   Import checkpoint size in number of variants\n"
	"  -C FLOAT Import checkpoint size in bases\n"
	"  -s       Hide all program messages [null]\n";
}

int import(int argc, char** argv){
	if(argc <= 2){
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

		{"checkpoint-variants",		optional_argument, 0,  'c' },
		{"checkpoint-bases",		optional_argument, 0,  'C' },
		{"permute",		no_argument, 0,  'p' },
		{"no-permute",		no_argument, 0,  'P' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	SILENT = 0;
	S32 checkpoint_n_variants = 500;
	double checkpoint_bp_window = 30e3;
	bool permute = true;

	while ((c = getopt_long(argc, argv, "i:o:c:C:spP?", long_options, &option_index)) != -1){
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
		case 'c':
			checkpoint_n_variants = atoi(optarg);
			if(checkpoint_n_variants <= 0){
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'C':
			checkpoint_bp_window = atof(optarg);
			if(checkpoint_bp_window <= 0){
				std::cerr << tachyon::utility::timestamp("ERROR") << "Cannot set checkpoint to <= 0..." << std::endl;
				return(1);
			}
			break;
		case 'p': permute = true;  break;
		case 'P': permute = false; break;
		case 's':
			SILENT = 1;
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
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling import..." << std::endl;
	}

	tachyon::VariantImporter importer(input, output, checkpoint_n_variants, checkpoint_bp_window);
	importer.setPermute(permute);

	if(!importer.Build())
		return 1;
	//importer.BuildCompressionDictionaries();

	return 0;
}
