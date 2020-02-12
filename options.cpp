#include <algorithm>
#include <unordered_set>
#include <iostream> 

#include <string.h>
#include <mpi.h>
#include <getopt.h>
#include <zlib.h>

#include "gottcha_db.h"
#include "options.h"
#include "mpi_util.h"
#include "deque_set.h"

using namespace std;

#define		KiloByte		size_t(1024)
#define		MegaByte		(KiloByte*KiloByte)
#define		GigaByte		(MegaByte*KiloByte)

#define 	DEFAULT_WORD_SIZE	24
#define		DEFAULT_MIN_FRAG	30

// A worker RAM value of 0 is a special value that directs each
// compute node to measure the amount of available RAM prior to
// initiating the signature calculations and use at most 
// RAM_FRACTION of this amount.
#define		DEFAULT_WORKER_RAM	0

// Make sure that we have been compiled against a version of zlib that
// supports writing uncompressed files
#if (ZLIB_VERNUM < 0x1280)
	#error "Please compile with a zlib version >= 1.2.8"
#endif

deque<string> split(const string &m_path, const char &m_delim);
size_t parse_memory_size(const string m_str);

Options::Options(int argc, char* argv[], bool m_load_options)
{	
	// Don't quit unless we need to
	quit = false;
	word_size = DEFAULT_WORD_SIZE;
	min_output_frag = DEFAULT_MIN_FRAG;
	verbose = Options::SILENT;
	compress_output = false;
	default_topology = UNKNOWN_TOPOLOGY;
	max_worker_RAM = DEFAULT_WORKER_RAM;
	
	if(m_load_options){

		// Command line options:
		// Taxa Level input and output
		// 	--strain <strain-level genome mapping file>
		// 	--strain.prefix <strain-level output file prefix>
		// 	--species <species-level genome mapping file>
		// 	--species.prefix <species-level output file prefix>
		// 	--genus <genus-level genome mapping file>
		// 	--genus.prefix <genus-level output file prefix>
		// 	--family <family-level genome mapping file>
		// 	--family.prefix <family-level output file prefix>
		// 	--order <order-level genome mapping file>
		// 	--order.prefix <order-level output file prefix>
		// 	--class <class-level genome mapping file>
		// 	--class.prefix <class-level output file prefix>
		// 	--phylum <phylum-level genome mapping file>
		// 	--phylum.prefix <phylum-level output file prefix>
		// 	--kingdom <kingdom-level genome mapping file>
		// 	--kingdom.prefix <kingdom-level output file prefix>
		// [--RAM <max RAM per MPI rank>]
		// [--squash <taxa label>]
		// [-w <digestion word size>]
		// [-f <min output fragment length>]
		// [--silent|--informative|--verbose]
		// [--compress (compress the output fasta files)]
		// [--linear (force all genomes to have linear topology)]
		// [--circular (force all genomes to have circular topology)]
		// [--log <log file>]
		// [--root <root directory for output files>]
		
		const char* options = "w:f:?h";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"squash", true, &config_opt, 1},
			{"silent", false, &config_opt, 2},
			{"informative", false, &config_opt, 3},
			{"verbose", false, &config_opt, 4},
			{"compress", false, &config_opt, 5},
			{"linear", false, &config_opt, 6},
			{"circular", false, &config_opt, 7},
			{"log", true, &config_opt, 8},
			{"root", true, &config_opt, 9},
			{"strain", true, &config_opt, 10},
			{"strain.prefix", true, &config_opt, 11},
			{"species", true, &config_opt, 12},
			{"species.prefix", true, &config_opt, 13},
			{"genus", true, &config_opt, 14},
			{"genus.prefix", true, &config_opt, 15},
			{"family", true, &config_opt, 16},
			{"family.prefix", true, &config_opt, 17},
			{"order", true, &config_opt, 18},
			{"order.prefix", true, &config_opt, 19},
			{"class", true, &config_opt, 20},
			{"class.prefix", true, &config_opt, 21},
			{"phylum", true, &config_opt, 22},
			{"phylum.prefix", true, &config_opt, 23},
			{"kingdom", true, &config_opt, 24},
			{"kingdom.prefix", true, &config_opt, 25},
			{"RAM", true, &config_opt, 26},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;
		
		bool print_usage = (argc == 1);
		
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){
	
			switch( opt_code ){
				case 0:
					if(config_opt == 1){ // squash
						
						groups_to_squash.push_back(optarg);
						break;
					}
						
					if(config_opt == 2){ // silent
					
						if(verbose < Options::SILENT){
							verbose = Options::SILENT;
						}

						break;
					}
					
					if(config_opt == 3){ // informative
					
						if(verbose < Options::INFORMATIVE){
							verbose = Options::INFORMATIVE;
						}
						
						break;
					}
					
					if(config_opt == 4){ // verbose
					
						if(verbose < Options::VERBOSE){
							verbose = Options::VERBOSE;
						}
						
						break;
					}
					
					if(config_opt == 5){ // compress
					
						compress_output = true;
						break;
					}
					
					if(config_opt == 6){ // linear
					
						default_topology = LINEAR;
						break;
					}
					
					if(config_opt == 7){ // circular
					
						default_topology = CIRCULAR;
						break;
					}
					
					if(config_opt == 8){ // log
					
						log_file = optarg;
						break;
					}
					
					if(config_opt == 9){ // root
					
						output_root_dir = optarg;
						break;
					}
					
					if(config_opt == 10){ // strain
					
						taxa_level_mapping_file[STRAIN_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 11){ // strain.prefix
					
						taxa_level_output_prefix[STRAIN_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 12){ // species
					
						taxa_level_mapping_file[SPECIES_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 13){ // species.prefix
					
						taxa_level_output_prefix[SPECIES_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 14){ // genus
					
						taxa_level_mapping_file[GENUS_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 15){ // genus.prefix
					
						taxa_level_output_prefix[GENUS_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 16){ // family
					
						taxa_level_mapping_file[FAMILY_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 17){ // family.prefix
					
						taxa_level_output_prefix[FAMILY_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 18){ // order
					
						taxa_level_mapping_file[ORDER_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 19){ // order.prefix
					
						taxa_level_output_prefix[ORDER_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 20){ // class
					
						taxa_level_mapping_file[CLASS_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 21){ // class.prefix
					
						taxa_level_output_prefix[CLASS_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 22){ // phylum
					
						taxa_level_mapping_file[PHYLUM_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 23){ // phylum.prefix
					
						taxa_level_output_prefix[PHYLUM_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 24){ // kingdom
					
						taxa_level_mapping_file[KINGDOM_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 25){ // kingdom.prefix
					
						taxa_level_output_prefix[KINGDOM_LEVEL] = optarg;
						break;
					}
					
					if(config_opt == 26){ // RAM.worker
					
						max_worker_RAM = parse_memory_size(optarg);
						break;
					}
					
					cerr << "Unknown flag!" << endl;
					break;
				case 'w':
					word_size = abs( atoi(optarg) );
					break;
				case 'f':
					min_output_frag = abs( atoi(optarg) );
					break;
				case 'h':
				case '?':
					print_usage = true;
					break;
				default:
					cerr << '\"' << (char)opt_code << "\" is not a valid option!" << endl;
					break;
			};
		}
		
		if(print_usage){
			
			quit = true;
			
			cerr << "Usage for gottcha_db (v. " << VERSION << "):" << endl;
			cerr << "\tTaxa Level input files and output prefixes" << endl;
			cerr << "\t\t--strain <strain-level genome mapping file>"<< endl;
			cerr << "\t\t--strain.prefix <strain-level output file prefix>"<< endl;
			cerr << "\t\t--species <species-level genome mapping file>"<< endl;
			cerr << "\t\t--species.prefix <species-level output file prefix>"<< endl;
			cerr << "\t\t--genus <genus-level genome mapping file>"<< endl;
			cerr << "\t\t--genus.prefix <genus-level output file prefix>"<< endl;
			cerr << "\t\t--family <family-level genome mapping file>"<< endl;
			cerr << "\t\t--family.prefix <family-level output file prefix>"<< endl;
			cerr << "\t\t--order <order-level genome mapping file>"<< endl;
			cerr << "\t\t--order.prefix <order-level output file prefix>"<< endl;
			cerr << "\t\t--class <class-level genome mapping file>"<< endl;
			cerr << "\t\t--class.prefix <class-level output file prefix>"<< endl;
			cerr << "\t\t--phylum <phylum-level genome mapping file>"<< endl;
			cerr << "\t\t--phylum.prefix <phylum-level output file prefix>"<< endl;
			cerr << "\t\t--kingdom <kingdom-level genome mapping file>"<< endl;
			cerr << "\t\t--kingdom.prefix <kingdom-level output file prefix>"<< endl;
			cerr << "\t[--RAM <max RAM per MPI rank>] (default is " << RAM_FRACTION*100 
				<< "% of available)" << endl;
			cerr << "\t[--squash <taxa label>] (prevent writing; can appear multiple times)" << endl;
			cerr << "\t[-w <digestion word size in bp>] (default is " << DEFAULT_WORD_SIZE << ")" << endl;
			cerr << "\t[-f <min output fragment length in bp>] (default is " << DEFAULT_MIN_FRAG << ")" << endl;
			cerr << "\t[--compress (compress the output fasta files)]"<< endl;
			cerr << "\t[--linear (force all genomes to have linear topology)]" << endl;
			cerr << "\t[--circular (force all genomes to have circular topology)]" << endl;
			cerr << "\t[--silent|--informative|--verbose] (default is slient)" << endl;
			cerr << "\t[--log <log file>]" << endl;
			cerr << "\t[--root <root directory for output files>]" << endl;
		}
		else{
			
			if( (word_size < 1) || (word_size > MAX_WORD_SIZE) ){

				cerr << "Please enter a valid word size (-w); 1 <= word size <= " << MAX_WORD_SIZE << endl;
				quit = true;
			}

			if(min_output_frag < 1){

				cerr << "Please enter a valid output minimum fragment size (-f); must be >= 1" << endl;
				quit = true;
			}

			// Make sure that the groups for which we will not be writing output for
			// appear once and only once
			make_set(groups_to_squash);
			
			// Make sure that *both* genome mapping files and prefixs have been provided
			// for all of the requested taxa-levels
			deque<TaxaLevel> requested_taxa;
			
			for(unordered_map<TaxaLevel, string>::const_iterator i = taxa_level_mapping_file.begin();
				i != taxa_level_mapping_file.end();++i){
				
				requested_taxa.push_back(i->first);
			}
			
			for(unordered_map<TaxaLevel, string>::const_iterator i = taxa_level_output_prefix.begin();
				i != taxa_level_output_prefix.end();++i){
				
				requested_taxa.push_back(i->first);
			}
			
			make_set(requested_taxa);
			
			for(deque<TaxaLevel>::const_iterator i = requested_taxa.begin();
				i != requested_taxa.end();++i){
				
				unordered_map<TaxaLevel, string>::const_iterator iter = 
					taxa_level_mapping_file.find(*i);
					
				if( iter == taxa_level_mapping_file.end() ){
				
					cerr << "Unable to find a " << taxa_level_name(*i) 
						<< "-level genome mapping file" << endl;
					quit = true;
				}
				
				iter = taxa_level_output_prefix.find(*i);
					
				if( iter == taxa_level_output_prefix.end() ){
				
					cerr << "Unable to find a " << taxa_level_name(*i) 
						<< "-level prefix" << endl;
					quit = true;
				}
			}

			if( requested_taxa.empty() ){
			
				cerr << "Please specify at least one taxa level input file and output prefix" << endl;
				quit = true;
			}
		}
	}
}

string taxa_level_name(const TaxaLevel &m_taxa)
{
	switch(m_taxa){
		case STRAIN_LEVEL:
			return "strain";
		case SPECIES_LEVEL:
			return "species";
		case GENUS_LEVEL:
			return "genus";
		case FAMILY_LEVEL:
			return "family";
		case ORDER_LEVEL:
			return "order";
		case CLASS_LEVEL:
			return "class";
		case PHYLUM_LEVEL:
			return "phylum";
		case KINGDOM_LEVEL:
			return "kingdom";
		default:
			throw __FILE__ ":taxa_level_name: Unknown taxa level";
	};
	
	throw __FILE__ ":taxa_level_name: Unknown taxa level";
	return string();
}

// Parse memory size as: \d+\w*[GB|MB|KB]
// and return the memory size in bytes
size_t parse_memory_size(const string m_str)
{
	// Extract the integer component
	string number;
	
	string::const_iterator iter = m_str.begin();
	
	while( ( iter != m_str.end() ) && isdigit(*iter) ){
		
		number.push_back(*iter);
		++iter;
	}
	
	// Skip any white space
	while( ( iter != m_str.end() ) && isspace(*iter) ){
		++iter;
	}
	
	// Convert the number to size_t
	double value = atof( number.c_str() );
	
	if(value < 0.0){
		throw __FILE__ ":parse_memory_size: Invalid input string!";
	}
	
	// Is there at least one more character?
	if( iter == m_str.end() ){
	
		// No, the number is in bytes by default
		return size_t(value);
	}
	
	switch( toupper(*iter) ){
		case 'B': // Bytes
			return size_t(value);
		case 'K': // Kilo Bytes
			return size_t(value*KiloByte);
		case 'M': // Mega Bytes
			return size_t(value*MegaByte);
		case 'G': // Giga Bytes
			return size_t(value*GigaByte);
	};
	
	// Should never get here
	throw __FILE__ ":parse_memory_size: Invalid input character detected";
	return 0;
}
