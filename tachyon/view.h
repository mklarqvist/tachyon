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
==============================================================================*/
#include <iostream>
#include <getopt.h>

#include <regex>

#include "support/enums.h"
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
	"  -O STRING output format: can be either JSON,VCF,BCF, or CUSTOM (-c must be triggered)\n"
	"  -f STRING interpreted filter string for slicing output (see manual)\n"
	"  -r STRING interval string\n"
	"  -R STRING path to file with interval strings\n"
	"  -d CHAR   output delimiter (-c must be triggered)\n"
	"  -y        custom output format (ignores VCF/BCF specification rules)\n"
	"  -G        drop all FORMAT fields from output\n"
	"  -h/H      header only / no header\n"
	"  -s        Hide all program messages\n\n"

	"Subset options:\n"
	"  -a, --trim-alt-alleles        trim alternate alleles not seen in the subset\n"
	"  -I, --no-update               do not (re)calculate INFO fields for the subset\n"
	"  -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n"
	"  -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n"
	"      --force-samples           only warn about unknown subset samples\n\n"

	"Filter options:\n"
    "  -c/C, --min-ac/--max-ac <int>[:<type>]      minimum/maximum count for non-reference (nref), 1st alternate (alt1), least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    "  -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes\n"
    "  -z/Z, --known/--novel                       select known/novel sites only (ID is not/is '.')\n"
    "  -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT (e.g. -m2 -M2 for biallelic sites)\n"
    "  -p/P, --phased/--exclude-phased             select/exclude sites where all samples are phased\n"
	"  -j,   --mixed-phasing                       select sites with both phased and unphased samples\n"
	"  -q/Q, --min-af/--max-af <float>[:<type>]    minimum/maximum frequency for non-reference (nref), 1st alternate (alt1), least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    "  -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype\n"
    "  -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]\n"
    "  -x/X, --private/--exclude-private           select/exclude sites where the non-reference alleles are exclusive (private) to the subset samples\n";
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
		{"output-type", optional_argument, 0,  'O' },
		{"vector-output", no_argument, 0,  'V' },
		{"annotate-genotype", no_argument, 0,  'X' },
		{"region",   optional_argument, 0,  'R' },
		{"noHeader",    no_argument, 0,  'H' },
		{"onlyHeader",  no_argument, 0,  'h' },
		{"dropFormat",  no_argument, 0,  'G' },
		{"customFormat",no_argument, 0,  'y' },
		{"silent",      no_argument, 0,  's' },
		{"af-min",      optional_argument, 0,  'q' },
		{"af-max",      optional_argument, 0,  'Q' },
		{"ac-min",      optional_argument, 0,  'c' },
		{"ac-max",      optional_argument, 0,  'C' },
		{"alleles-min",      optional_argument, 0,  'm' },
		{"alleles-max",      optional_argument, 0,  'M' },
		{"known",      no_argument, 0,  'z' },
		{"novel",      no_argument, 0,  'Z' },
		{"phased",      no_argument, 0,  'p' },
		{"exclude-phased",      no_argument, 0,  'P' },
		{"mixed-phase",      no_argument, 0,  'j' },
		{"uncalled",      no_argument, 0,  'u' },
		{"exclude-uncalled",      no_argument, 0,  'U' },
		{0,0,0,0}
	};

	tachyon::VariantReaderSettings settings;
	tachyon::DataBlockSettings block_settings;
	tachyon::VariantReaderFilters filters;
	std::vector<std::string> interpret_commands;

	SILENT = 0;
	std::string temp;

	while ((c = getopt_long(argc, argv, "i:o:k:f:d:O:R:yGshHVX?q:Q:m:M:pPuUc:C:jzZ", long_options, &option_index)) != -1){
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
		case 'k':
			settings.keychain_file = std::string(optarg);
			break;
		case 'q':
			filters.filter_af(atof(optarg));
			break;
		case 'Q':
			filters.filter_af(atof(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;

		case 'm':
			filters.filter_n_alts(atoi(optarg));
			break;
		case 'M':
			filters.filter_n_alts(atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;

		case 'c':
			filters.filter_ac(atoi(optarg));
			break;
		case 'C':
			filters.filter_ac(atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;

		case 'p':
			filters.filter_uniform_phase(true);
			break;
		case 'P':
			filters.filter_uniform_phase(false);
			break;
		case 'j':
			filters.filter_mixed_phase(true);
			break;

		case 'u':
			filters.filter_missing(true);
			break;
		case 'U':
			filters.filter_missing(false);
			break;

		case 'z':
			filters.filter_known_novel(true);
			break;
		case 'Z':
			filters.filter_known_novel(false);
			break;

		case 'f':
			interpret_commands.push_back(std::string(optarg));
			break;
		case 'G':
			settings.drop_format = true;
			break;
		case 's':
			SILENT = 1;
			break;
		case 'y':
			settings.custom_output_format = true;
			break;
		case 'R':
			settings.interval_strings.push_back(std::string(optarg));
			break;
		case 'd':
			settings.custom_delimiter = true;
			temp = std::string(optarg);
			if(temp.size() != 1 && !(temp[0] == '\\' && temp.size() == 2)){
				std::cerr << "not a legal delimiter" << std::endl;
				return(1);
			}
			if(temp.size() == 1) settings.custom_delimiter_char = temp[0];
			else {
				if(temp[1] == 't') settings.custom_delimiter_char = '\t';
				else if(temp[1] == 'n') settings.custom_delimiter_char = '\n';
				else {
					std::cerr << "not a legal delimiter" << std::endl;
					return(1);
				}
			}
			break;
		case 'h':
			settings.header_only = true;
			break;
		case 'H':
			settings.show_header = false;
			break;
		case 'V':
			settings.output_FORMAT_as_vector = true;
			break;

		case 'O':
			settings.output_type = std::string(optarg);
			break;

		case 'X':
			settings.annotate_genotypes = true;
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
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling view..." << std::endl;
	}

	tachyon::VariantReader reader;
	reader.getSettings() = settings;

	// temp
	if(settings.keychain_file.size()){
		std::ifstream keychain_reader(settings.keychain_file, std::ios::binary | std::ios::in);
		if(!keychain_reader.good()){
			std::cerr << tachyon::utility::timestamp("ERROR") <<  "Failed to open keychain: " << settings.keychain_file << "..." << std::endl;
			return 1;
		}

		keychain_reader >> reader.keychain;
		if(!keychain_reader.good()){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse keychain..." << std::endl;
			return 1;
		}
	}

	if(!reader.open(settings.input)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << settings.input << "..." << std::endl;
		return 1;
	}

	if(settings.header_only){
		reader.header.literals += "\n##tachyon_viewVersion=" + tachyon::constants::PROGRAM_NAME + "-" + VERSION + ";";
		reader.header.literals += "libraries=" +  tachyon::constants::PROGRAM_NAME + '-' + tachyon::constants::TACHYON_LIB_VERSION + ","
				  + SSLeay_version(SSLEAY_VERSION) + "," + "ZSTD-" + ZSTD_versionString() + "; timestamp=" + tachyon::utility::datetime();

		reader.header.literals += "\n##tachyon_viewCommand=" + tachyon::constants::LITERAL_COMMAND_LINE + "\n";
		reader.header.literals += reader.getSettings().get_settings_string();

		std::cout << reader.header.literals << std::endl;
		reader.header.writeVCFHeaderString(std::cout, true);
		return(0);
	}

	// User provided '-f' string(s)
	if(interpret_commands.size()){
		if(!reader.getBlockSettings().parseCommandString(interpret_commands, reader.header, settings.custom_output_format)){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse command..." << std::endl;
			return(1);
		}
	} else {
		reader.getBlockSettings().loadAll(true);

		if(settings.drop_format){
			reader.getBlockSettings().loadGenotypes(false);
			reader.getBlockSettings().load_ppa    = false;
			reader.getBlockSettings().load_format = false;
		}
	}

	if(settings.custom_delimiter){
		if(settings.custom_output_format == false){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Have to trigger -c when using a custom separator" << std::endl;
			return(1);
		}
		reader.getBlockSettings().setCustomDelimiter(settings.custom_delimiter_char);
	}

	reader.getBlockSettings().output_format_vector = settings.output_FORMAT_as_vector;

	if(settings.output_type.size()){
		std::transform(settings.output_type.begin(), settings.output_type.end(), settings.output_type.begin(), ::toupper); // transform to UPPERCASE
		if(strncmp(&settings.output_type[0], "JSON", 4) == 0 && settings.output_type.size() == 4){
			if(settings.custom_delimiter)
				std::cerr << tachyon::utility::timestamp("WARNING") << "Custom output delimiter is incompatible with JSON. Disabled..." << std::endl;

			settings.custom_output_format = true;
			reader.getBlockSettings().custom_output_format = true;
			reader.getBlockSettings().output_json = true;
			reader.getBlockSettings().output_format_vector = true;
		} else if(strncmp(&settings.output_type[0], "VCF", 3) == 0 && settings.output_type.size() == 3){
			if(settings.custom_delimiter)
				std::cerr << tachyon::utility::timestamp("WARNING") << "Custom output delimiter is incompatible with VCF. Disabled..." << std::endl;

			reader.getBlockSettings().custom_output_format = false;
			reader.getBlockSettings().custom_delimiter = false;
			reader.getBlockSettings().custom_delimiter_char = '\t';

			if(settings.output_FORMAT_as_vector)
				std::cerr << tachyon::utility::timestamp("WARNING") << "Output FORMAT as vectors (-V) is incompatible with VCF output. Disabled..." << std::endl;

			reader.getBlockSettings().output_format_vector = false;
		} else if(strncmp(&settings.output_type[0], "BCF", 3) == 0 && settings.output_type.size() == 3){
			reader.getBlockSettings().custom_output_format = false;
			std::cerr << tachyon::utility::timestamp("ERROR") << "BCF output not supported yet." << std::endl;
			return(1);
		} else if(strncmp(&settings.output_type[0], "CUSTOM", 6) == 0 && settings.output_type.size() == 6){
			reader.getBlockSettings().custom_output_format = true;
			settings.custom_output_format = true;
		} else {
			std::cerr << tachyon::utility::timestamp("ERROR") << "Unrecognised output option: " << settings.output_type << "..." << std::endl;
			return(1);
		}
	}

	// If user is triggering annotation
	if(settings.annotate_genotypes){
		reader.getBlockSettings().annotate_extra = true;
		reader.getBlockSettings().loadGenotypes(true);
		reader.getBlockSettings().load_set_membership = true;
		reader.getBlockSettings().load_alleles = true;
		reader.getBlockSettings().load_positons = true;;
	}

	reader.getFilterSettings() = filters;

	tachyon::algorithm::Timer timer;
	timer.Start();

	if(settings.show_header) reader.getBlockSettings().show_vcf_header = true;
	else reader.getBlockSettings().show_vcf_header = false;

	// Todo: pass to settings

	//reader.getSettings().interval_strings.push_back("20");
	//reader.getSettings().interval_strings.push_back(" 20:1e6   ");
	//reader.getSettings().interval_strings.push_back(" 20: 1e6 -2e6");
	if(reader.getSettings().validateIntervalStrings() == false) return(1);
	//return(1);
	if(reader.parseIntervals() == false) return(1);

	// Temp
	//while(reader.nextBlock()) reader.getGenotypeSummary(std::cout);
	//return(0);

	U64 n_variants = 0;
	if(settings.custom_output_format) n_variants = reader.outputCustom();
	else n_variants = reader.outputVCF();

	//std::cerr << "Blocks: " << n_blocks << std::endl;
	std::cerr << "Variants: " << tachyon::utility::ToPrettyString(n_variants) << " genotypes: " << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t' << timer.ElapsedString() << '\t' << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count())) << std::endl;

	return 0;
}
