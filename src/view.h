/*
Copyright (C) 2017-current Genome Research Ltd.
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
	"About:  Convert YON->VCF/BCF/\n"
	"Usage:  " << Tachyon::Constants::PROGRAM_NAME << " view [options] -i <in.vcf>/<in.bcf> -o <output>\n\n"
	"Options:\n"
	"  -i FILE  input Tachyon file (required)\n"
	"  -o FILE  output file prefix or - for stdout\n"
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
		{"dropFormat",	no_argument, 0,  'G' },
		{"dropInfo",	no_argument, 0,  'I' },
		{"silent",		no_argument, 0,  's' },
		{0,0,0,0}
	};

	std::string input;
	std::string output;
	SILENT = 0;
	bool dropFormat = false;
	bool dropInfo   = false;

	while ((c = getopt_long(argc, argv, "i:o:GIs?", long_options, &option_index)) != -1){
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
		case 'G':
			dropFormat = true;
			break;
		case 'I':
			dropInfo = true;
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

	Tachyon::TachyonReader reader;

	reader.settings.loadAll();
	//reader.settings.loadPPA = false;
	if(dropInfo)   reader.settings.unsetInfo();
	if(dropFormat) reader.settings.unsetFormat();

	//reader.settings.loadMetaHot = true;
	//reader.settings.loadGT = true;
	//reader.settings.loadGTSimple = true;


	//reader.settings.loadMetaHot = true;
	//reader.settings.loadMetaCold = true;
	// Todo: deduplicate and move to function in settings class
	//for(U32 i = 14; i < 18; ++i)
	//	reader.settings.load_info_ID.push_back(i);
	//
	/*

	reader.settings.load_info_ID.push_back(5);
	reader.settings.load_info_ID.push_back(15);
	reader.settings.load_info_ID.push_back(16);
	reader.settings.load_info_ID.push_back(17);
	reader.settings.load_info_ID.push_back(18);
	reader.settings.load_info_ID.push_back(11);
	reader.settings.load_info_ID.push_back(24);
	reader.settings.load_info_ID.push_back(11);
	reader.settings.load_info_ID.push_back(25);
	*/

	if(!reader.open(input)){
		std::cerr << "failed to open" << std::endl;
		return 1;
	}

	while(reader.getNextBlock()){
		//reader.toVCFStringFast();
		//reader.toVCFString();
		//reader.iterateMeta();
		reader.iterateGT();
	}

	return 0;
}
