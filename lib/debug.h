#ifndef DEBUG_H_
#define DEBUG_H_

#include <iostream>
#include <getopt.h>

#include "utility.h"
#include "variant_reader.h"

void debug_usage(void){

}

int debug(int argc, char** argv){
	if(argc <= 2){
		programHelp();
		return(1);
	}

	int c;
	if(argc < 2){
		debug_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",    required_argument, 0, 'i' },
		{"output",   optional_argument, 0, 'o' },
		{"threads",  optional_argument, 0, 't' },
		{"silent",   no_argument,       0, 's' },
		{0,0,0,0}
	};

	tachyon::VariantReaderSettings settings;
	SILENT = 0;
	int n_threads = 1;

	while ((c = getopt_long(argc, argv, "i:o:t:s?", long_options, &option_index)) != -1){
		switch (c){
		case 0:
			std::cerr << "Case 0: " << option_index << '\t' << long_options[option_index].name << std::endl;
			break;
		case 'i':
			settings.input = std::string(optarg);
			break;
		case 'o':
			settings.output = std::string(optarg);
			break;
		case 't':
			n_threads = atoi(optarg);
			if(n_threads <= 0) return 1;
			break;
		case 's':
			SILENT = 1;
			break;
		default:
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognized option: " << (char)c << std::endl;
			return(1);
		}
	}

	if(settings.input.length() == 0){
		std::cerr << tachyon::utility::timestamp("ERROR") << "No input value specified..." << std::endl;
		return(1);
	}

	// Print messages
	if(!SILENT){
		programMessage();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling debug..." << std::endl;
	}

	tachyon::VariantReader reader;
	reader.GetSettings() = settings;
	//reader.GetBlockSettings().LoadAll(true).LoadAllFormat(false);
	reader.GetBlockSettings().LoadAll(true);

	if(!reader.open(settings.input)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << settings.input << "..." << std::endl;
		return 1;
	}

	reader.TempWrite();

	return 0;
}



#endif /* DEBUG_H_ */
