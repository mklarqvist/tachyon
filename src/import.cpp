/*
Copyright (C) 2016-2017 Genome Research Ltd.
Author: Marcus D. R. Klarqvist <mk21@sanger.ac.uk>

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
#include "core/Importer.h"

void import_usage(void){
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

int main(int argc, char** argv){
	if(Tachyon::Helpers::isBigEndian()){
		std::cerr << Tachyon::Helpers::timestamp("ERROR") << "Tomahawk does not support big endian systems..." << std::endl;
		return(1);
	}

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
		{"checkpoint",		optional_argument, 0,  'c' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	SILENT = 0;
	S32 checkpoint = 500;

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
		case 'c':
			checkpoint = atoi(optarg);
			if(checkpoint < 0){
				std::cerr << Tachyon::Helpers::timestamp("ERROR") << "Cannot set checkpoint to < 0..." << std::endl;
				return(1);
			}
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
		std::cerr << Tachyon::Helpers::timestamp("LOG") << "Calling import..." << std::endl;
	}

	//Tomahawk::TomahawkImportWriter writer;
	//writer.test();
	/*
	Tomahawk::Index::IndexBlockEntryBase b;
	std::cerr << b.contigID << std::endl;
	Tomahawk::Index::IndexBlockEntry e;
	std::cerr << e.contigID << std::endl;
	return(1);
	*/


	Tachyon::Importer importer(input, output, checkpoint);

	if(!importer.Build())
		return 1;

	return 0;
}
