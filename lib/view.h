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
	"  -y        custom output format (ignores VCF specification rules)\n"
	"  -V        custom output data as vectors instead of per sample (valid only with -y)\n"
	"  -G        drop all FORMAT fields from output\n"
	"  -h/H      header only / no header\n"
	"  -s        Hide all program messages\n\n"

	"Subset options:\n"
	"  -s, --samples [^]<list>       comma separated list of samples to include (or exclude with \"^\" prefix)\n"
	"  -S, --samples-file [^]<file>  file of samples to include (or exclude with \"^\" prefix)\n\n"

	"Filter options:\n"
	"  -a/A, --ref-match/--alt-match <REGEX>       regular expression pattern for the reference allele -a or for any alternative alleles -A\n"
	"  -n,   --name-match <REGEX>                  regular expression pattern for the locus name\n"
    "  -c/C, --min-ac/--max-ac <int>               minimum/maximum count for non-reference least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    "  -g,   --genotype [^]<hom|het|miss>          require one or more hom/het/missing genotype or, if prefixed with \"^\", exclude sites with hom/het/missing genotypes\n"
    "  -z/Z, --known/--novel                       select known/novel sites only (ID is not/is '.')\n"
	"  -q/Q, --min-quality/--max-quality           minimum/maximum quality value\n"
    "  -m/M, --min-alleles/--max-alleles <int>     minimum/maximum number of alleles listed in REF and ALT\n"
    "  -p/P, --phased/--exclude-phased             select/exclude sites where all samples are phased\n"
	"  -j/J  --mixed-phasing/--no-mixed-phasing    select sites with both phased and unphased samples\n"
	"  -w/W  --mixed-ploidy/--no-mixed-ploidy      select sites with mixed ploidy\n"
	"  -l/L, --min-af/--max-af <float>             minimum/maximum frequency for non-reference least frequent\n"
    "                                                 (minor), most frequent (major) or sum of all but most frequent (nonmajor) alleles [nref]\n"
    "  -u/U, --uncalled/--exclude-uncalled         select/exclude sites without a called genotype\n"
	"  -e/E, --remove-unseen/--keep-unseen         select/exclude sites with unseen alternative allele(s)\n"
    "  -v/V, --types/--exclude-types <list>        select/exclude comma-separated list of variant types: snps,indels,mnps,ref,bnd,other [null]\n\n";
}

int view(int argc, char** argv){
	if(argc < 2){
		programHelp();
		return(1);
	}

	int c;
	if(argc == 2){
		view_usage();
		return(1);
	}

	int option_index = 0;
	static struct option long_options[] = {
		{"input",             required_argument, 0,  'i' },
		{"output",            optional_argument, 0,  'o' },
		{"keychain",          optional_argument, 0,  'k' },
		{"filter",            optional_argument, 0,  'f' },
		{"delimiter",         optional_argument, 0,  'd' },
		{"output-type",       optional_argument, 0,  'O' },
		{"vector-output",     no_argument,       0,  'V' },
		{"annotate-genotype", no_argument,       0,  'X' },
		{"region",            optional_argument, 0,  'r' },
		{"noHeader",          no_argument,       0,  'H' },
		{"onlyHeader",        no_argument,       0,  'h' },
		{"dropFormat",        no_argument,       0,  'G' },
		{"customFormat",      no_argument,       0,  'y' },
		{"silent",            no_argument,       0,  's' },
		{"af-min",            optional_argument, 0,  'l' },
		{"af-max",            optional_argument, 0,  'L' },
		{"ac-min",            optional_argument, 0,  'c' },
		{"ac-max",            optional_argument, 0,  'C' },
		{"alleles-min",       optional_argument, 0,  'm' },
		{"alleles-max",       optional_argument, 0,  'M' },
		{"known",             no_argument,       0,  'z' },
		{"novel",             no_argument,       0,  'Z' },
		{"phased",            no_argument,       0,  'p' },
		{"exclude-phased"    ,no_argument,       0,  'P' },
		{"mixed-phase",       no_argument,       0,  'j' },
		{"no-mixed-phase",    no_argument,       0,  'J' },
		{"uncalled",          no_argument,       0,  'u' },
		{"exclude-uncalled",  no_argument,       0,  'U' },
		{"ref-match",         optional_argument, 0,  'a' },
		{"alt-match",         optional_argument, 0,  'A' },
		{"name-match",        optional_argument, 0,  'n' },
		{"mixed-ploidy",      no_argument,       0,  'w' },
		{"no-mixed-ploidy",   no_argument,       0,  'W' },
		{"remove-unseen",     no_argument,       0,  'e' },
		{"keep-unseen",       no_argument,       0,  'E' },
		{"min-quality",       optional_argument, 0,  'q' },
		{"max-quality",       optional_argument, 0,  'Q' },
		{0,0,0,0}
	};

	tachyon::VariantReaderSettings settings;
	tachyon::DataBlockSettings     block_settings;
	std::vector<std::string>       interpret_commands;
	std::vector<std::string>       interval_strings;

	SILENT = 0;
	std::string temp;
	tachyon::VariantReader reader;
	tachyon::VariantReaderFilters& filters = reader.getFilterSettings();

	while ((c = getopt_long(argc, argv, "i:o:k:f:d:O:r:yGshHVX?l:L:m:M:pPuUc:C:jJzZa:A:n:wWeEq:Q:", long_options, &option_index)) != -1){
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
		case 'l':
			filters.add(tachyon::YON_FILTER_ALLELE_FREQUENCY, atof(optarg), tachyon::YON_CMP_GREATER);
			filters.require_genotypes = true;
			break;
		case 'L':
			filters.add(tachyon::YON_FILTER_ALLELE_FREQUENCY, atof(optarg), tachyon::YON_CMP_LESS_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'q':
			filters.add(tachyon::YON_FILTER_QUALITY, atof(optarg), tachyon::YON_CMP_GREATER);
			break;
		case 'Q':
			filters.add(tachyon::YON_FILTER_QUALITY, atof(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'm':
			filters.add(tachyon::YON_FILTER_NUMBER_ALT_ALLELES, atoi(optarg), tachyon::YON_CMP_GREATER);
			break;
		case 'M':
			filters.add(tachyon::YON_FILTER_NUMBER_ALT_ALLELES, atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			break;
		case 'a':
			filters.add(tachyon::YON_FILTER_REFERENCE_ALLELE, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 'A':
			filters.add(tachyon::YON_FILTER_ALT_ALLELE, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 'n':
			filters.add(tachyon::YON_FILTER_NAME, std::string(optarg), tachyon::YON_CMP_REGEX);
			break;
		case 'c':
			filters.add(tachyon::YON_FILTER_ALLELE_COUNT, atoi(optarg), tachyon::YON_CMP_GREATER);;
			filters.require_genotypes = true;
			break;
		case 'C':
			filters.add(tachyon::YON_FILTER_ALLELE_COUNT, atoi(optarg), tachyon::YON_CMP_LESS_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'p':
			filters.add(tachyon::YON_FILTER_UNIFORM_PHASE, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'P':
			filters.add(tachyon::YON_FILTER_UNIFORM_PHASE, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'j':
			filters.add(tachyon::YON_FILTER_MIXED_PHASING, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'J':
			filters.add(tachyon::YON_FILTER_MIXED_PHASING, (bool)false, tachyon::YON_CMP_EQUAL);
			break;
		case 'w':
			filters.add(tachyon::YON_FILTER_MIXED_PLOIDY, (bool)true, tachyon::YON_CMP_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'W':
			filters.add(tachyon::YON_FILTER_MIXED_PLOIDY, (bool)false, tachyon::YON_CMP_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'u':
			filters.add(tachyon::YON_FILTER_MISSING_GT, 0, tachyon::YON_CMP_GREATER);
			filters.require_genotypes = true;
			break;
		case 'U':
			filters.add(tachyon::YON_FILTER_MISSING_GT, 0, tachyon::YON_CMP_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'z':
			filters.add(tachyon::YON_FILTER_KNOWN_NOVEL, (bool)true, tachyon::YON_CMP_EQUAL);
			break;
		case 'Z':
			filters.add(tachyon::YON_FILTER_KNOWN_NOVEL, (bool)false, tachyon::YON_CMP_EQUAL);
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
		case 'r':
			interval_strings.push_back(std::string(optarg));
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
		case 'e':
			filters.add(tachyon::YON_FILTER_UNSEEN_ALT, (bool)true, tachyon::YON_CMP_EQUAL);
			filters.require_genotypes = true;
			break;
		case 'E':
			filters.add(tachyon::YON_FILTER_UNSEEN_ALT, (bool)false, tachyon::YON_CMP_EQUAL);
			filters.require_genotypes = true;
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
	/*
	if(!SILENT){
		programHelp();
		std::cerr << tachyon::utility::timestamp("LOG") << "Calling view..." << std::endl;
	}
	*/


	reader.getSettings() = settings;

	if(!reader.open(settings.input)){
		std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to open file: " << settings.input << "..." << std::endl;
		return 1;
	}

	if(settings.header_only){
		reader.printHeaderVCF();
		return(0);
	}

	// User provided '-f' string(s)
	if(interpret_commands.size()){
		if(!reader.getBlockSettings().parseCommandString(interpret_commands, reader.getGlobalHeader(), settings.custom_output_format)){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Failed to parse command..." << std::endl;
			return(1);
		}
	} else {
		reader.getBlockSettings().loadAll(true);

		if(settings.drop_format){
			reader.getBlockSettings().loadGenotypes(false);
			reader.getBlockSettings().ppa(false, false);
			reader.getBlockSettings().format_all(false, false);
		}
	}

	if(settings.custom_delimiter){
		if(settings.custom_output_format == false){
			std::cerr << tachyon::utility::timestamp("ERROR") << "Have to trigger -y when using a custom separator" << std::endl;
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
		reader.getBlockSettings().set_membership(true, true);
		reader.getBlockSettings().alleles(true, true);
		reader.getBlockSettings().positions(true, true);
	}

	if(filters.doRequireGenotypes()){
		reader.getBlockSettings().loadGenotypes(true);
		reader.getBlockSettings().set_membership.load = true;
		reader.getBlockSettings().alleles.load = true;
		reader.getBlockSettings().positions.load = true;
	}

	reader.getBlockSettings().parseSettings(reader.getGlobalHeader());

	tachyon::algorithm::Timer timer;
	timer.Start();

	if(settings.show_header) reader.getBlockSettings().show_vcf_header = true;
	else reader.getBlockSettings().show_vcf_header = false;

	if(reader.addIntervals(interval_strings) == false) return(1);

	U64 n_variants = 0;
	if(settings.custom_output_format) n_variants = reader.outputCustom();
	else n_variants = reader.outputVCF();

	//std::cerr << "Blocks: " << n_blocks << std::endl;
	/*
	std::cerr << "Variants: "
	          << tachyon::utility::ToPrettyString(n_variants) << " genotypes: "
	          << tachyon::utility::ToPrettyString(n_variants*reader.header.getSampleNumber()) << '\t'
	          << timer.ElapsedString() << '\t'
	          << tachyon::utility::ToPrettyString((U64)((double)n_variants*reader.header.getSampleNumber()/timer.Elapsed().count()))
	          << std::endl;
	*/

	return 0;
}
