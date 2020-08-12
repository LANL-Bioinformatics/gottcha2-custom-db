// Create the input file for the gottcha_db program
// Inputs: 
//	1) A directory of GenBank flat files (gbk/gb/gbff) to read
//	2) The NCBI taxonomy files (nodes.dmp, names.dmp and gi_taxid_nucl.dmp)
//	3) The desired taxonomic level to output (kingdom, phylum, class, order, family, genus, species, strain)
// Outputs:
//	1) The tab-delimited file containing unique taxa_id strings, file locations
//	   and molecule topology (linear/circular) information.
//
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Thu Nov 19 09:14:06 2015

#include <iostream>
#include <fstream>
#include <sstream>
#include <deque>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <zlib.h>

#include "accession.h"
#include "sequence.h"

// File manipulation for linux (this follows the 4.2 and 4.3 BSD standard -- see
// "Using C on the Unix System" by D. Curry).
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

using namespace std;

#define		VERSION			"1.3"
#define		PATH_SEPERATOR		"/"

// GenBank flat files and fasta are currently supported (should be easy to
// support EMBL files too).
struct FileInfo {

	const char *ext;
	FileType file_type;	
};

const FileInfo file_types[] = {
	{".gbk", GENBANK}, {".gbk.gz", GENBANK},
	{".gbff", GENBANK}, {".gbff.gz", GENBANK},
	{".gb", GENBANK}, {".gb.gz", GENBANK},
	{".fa", FASTA}, {".fa.gz", FASTA},
	{".fna.gz", FASTA}, {".fna.gz", FASTA},
	{NULL, UNKNOWN_FILETYPE}
};

// Enumerate all of the valid NCBI taxonomic levels + a special level called "STRAIN"
typedef enum {SUPERKINGDOM, KINGDOM, SUBKINGDOM, SUPERPHYLUM, PHYLUM, SUBPHYLUM, 
	SUPERCLASS, CLASS, SUBCLASS, INFRACLASS, SUPERORDER, ORDER, 
	PARVORDER, SUBORDER, INFRAORDER, SUPERFAMILY, FAMILY, 
	SUBFAMILY, GENUS, SUBGENUS, TRIBE, SUBTRIBE,
	SPECIES, SPECIES_GROUP, SPECIES_SUBGROUP,
	SUBSPECIES, STRAIN, NO_RANK, 
	VARIETAS, FORMA, UNKNOWN} Taxonomy;

typedef int TaxID;

struct TreeNode
{
	TaxID child;
	TaxID parent;
	Taxonomy rank;
	
	TreeNode()
	{
		child = 0;
		parent = 0;
		rank = UNKNOWN;
	};

	// Only sort by the child taxid
	inline bool operator<(const TreeNode &m_rhs) const
	{
		return (child < m_rhs.child);
	};

	inline bool operator<(const TaxID &m_child) const
	{
		return (child < m_child);
	};	
};

struct SeqInfo
{
	SeqInfo()
	{
		accession = 0;
		topo = UNKNOWN_TOPOLOGY;
	};
	
	Accession accession;
	Topology topo;
};

struct find_by_accession
{
	inline bool operator()(const pair<Accession, TaxID> &m_a, const Accession &m_b) const
	{
		return (m_a.first < m_b);
	};
};

void find_sequence(deque<string> &m_data, const string &m_path);
void parse_names(unordered_map<TaxID, string> &m_names, const string &m_file);
void parse_nodes(deque<TreeNode> &m_tree, const string &m_file);
void parse_accession_to_taxid(deque< pair<Accession, TaxID> > &m_data, const string &m_file, 
	const deque<Accession> &m_target);
void parse_sequence(unordered_multimap<string, SeqInfo> &m_info, const string &m_file);
void parse_sequence_gbk(unordered_multimap<string, SeqInfo> &m_info, const string &m_file);
void parse_sequence_fasta(unordered_multimap<string, SeqInfo> &m_info, const string &m_file);
int string_to_int(const string &m_buffer);
string trim_white_space(const string &m_string);
string tolower(const string &m_str);
string remove_whitespace(const string &m_str);
deque<string> split(const string &m_path, const char &m_delim);
TaxID get_output_taxa(const TaxID &m_id, const deque<TreeNode> &m_tree, 
	const Taxonomy &m_level);
FileType get_file_type(const string &m_file);

int main(int argc, char *argv[])
{
	try{
		
		// Command line options
		// -o <tab delimited output file *prefix*>
		// --ncbi <ncbi parent directory with default file names>
		// [--kingdom | --phylum | --class | --order | --family | --genus | --species | --strain]
		// [--ncbi.names <read as names.dmp file>]
		// [--ncbi.nodes <read as nodes.dmp file>]
		// [--ncbi.accession <NCBI accession to taxid file> (may be repeated)]
		// [--ignore-missing-taxid] (ignore accessions with missing taxid information)
		// [--compress] (compress the output files with zlib/gzip)
		// <directories or files to search for sequence files> ... (may be repeated multiple times)
		const char* options = "o:?h";
		int config_opt = 0;
		int long_index = 0;

		struct option long_opts[] = {
			{"ncbi", true, &config_opt, 1},
			{"ncbi.names", true, &config_opt, 2},
			{"ncbi.nodes", true, &config_opt, 3},
			{"ncbi.accession", true, &config_opt, 4},
			{"kingdom", false, &config_opt, 5},
			{"phylum", false, &config_opt, 6},
			{"class", false, &config_opt, 7},
			{"order", false, &config_opt, 8},
			{"family", false, &config_opt, 9},
			{"genus", false, &config_opt, 10},
			{"species", false, &config_opt, 11},
			{"strain", false, &config_opt, 12},
			{"ignore-missing-taxid", false, &config_opt, 13},
			{"superkingdom", false, &config_opt, 14},
			{"compress", false, &config_opt, 15},
			{0,0,0,0} // Terminate options list
		};

		int opt_code;
		opterr = 0;
		
		bool print_usage = (argc == 1);
		string output_file_prefix;
		string ncbi_dir;
		string ncbi_names_file;
		string ncbi_nodes_file;
		deque<string> ncbi_accession2taxid_file;
		bool ignore_missing_taxid = false;
		deque<Taxonomy> t;
		bool compress = false;
		
		while( (opt_code = getopt_long( argc, argv, options, long_opts, &long_index) ) != EOF ){
	
			switch( opt_code ){
				case 0:
					if(config_opt == 1){ // ncbi
						
						ncbi_dir = optarg;
						break;
					}
						
					if(config_opt == 2){ // ncbi.names
					
						ncbi_names_file = optarg;
						break;
					}
					
					if(config_opt == 3){ // ncbi.nodes
					
						ncbi_nodes_file = optarg;
						break;
					}
					
					if(config_opt == 4){ // ncbi.accession
					
						ncbi_accession2taxid_file.push_back(optarg);
						break;
					}
					
					if(config_opt == 5){ // kingdom
					
						t.push_back(KINGDOM);
						break;
					}
					
					if(config_opt == 6){ // phylum
					
						t.push_back(PHYLUM);
						break;
					}
					
					if(config_opt == 7){ // class
					
						t.push_back(CLASS);
						break;
					}

					if(config_opt == 8){ // order
					
						t.push_back(ORDER);
						break;
					}
					
					if(config_opt == 9){ // family
					
						t.push_back(FAMILY);
						break;
					}
					
					if(config_opt == 10){ // genus
					
						t.push_back(GENUS);
						break;
					}
					
					if(config_opt == 11){ // species

						t.push_back(SPECIES);
						break;
					}
					
					if(config_opt == 12){ // strain
					
						t.push_back(STRAIN);
						break;
					}
					
					if(config_opt == 13){ // ignore-missing-taxid
					
						ignore_missing_taxid = true;
						break;
					}
					
					if(config_opt == 14){ // superkingdom
					
						t.push_back(SUPERKINGDOM);
						break;
					}
					
					if(config_opt == 15){ // compress
					
						compress = true;
						break;
					}
					
					cerr << "Unknown flag!" << endl;
					break;
				case 'o':
					output_file_prefix = optarg;
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
			
			cerr << "Usage for gottcha_format (v. " << VERSION << "):" << endl;
			cerr << "\t-o <output file prefix>" << endl;
			cerr << "\t--ncbi <ncbi parent directory with default file names>" << endl;
			cerr << "\tTaxonomic target (select any combination of):" << endl;
			cerr << "\t\t--superkingdom" << endl;
			cerr << "\t\t--kingdom" << endl;
			cerr << "\t\t--phylum" << endl;
			cerr << "\t\t--class" << endl;
			cerr << "\t\t--order" << endl;
			cerr << "\t\t--family" << endl;
			cerr << "\t\t--genus" << endl;
			cerr << "\t\t--species" << endl;
			cerr << "\t\t--strain" << endl;
			cerr << "\t[--ncbi.names <read as names.dmp file>] (overrides default)"<< endl;
			cerr << "\t[--ncbi.nodes <read as nodes.dmp file>] (overrides default)" << endl;
			cerr << "\t[--ncbi.accession <NCBI accession to taxid file>]" << endl;
			cerr << "\t\tDefault file = <ncbi>/nucl_gb.accession2taxid.gz" << endl;
			cerr << "\t\tDefault file = <ncbi>/nucl_wgs.accession2taxid.gz" << endl;
			cerr << "\t[--ignore-missing-taxid] (ignore accessions with missing taxid information)" << endl;
			cerr << "\t[--compress] (compress the output files with zlib/gzip)" << endl;
			cerr << "\t<directory or file to search for sequence files> ... (may be repeated multiple times)" << endl;
			
			return EXIT_SUCCESS;
		}
		
		if( t.empty() ){
			
			cerr << "Please specify at least one taxonomic level to group sequences" << endl;
			return EXIT_FAILURE;
		}
				
		if( output_file_prefix.empty() ){
			
			cerr << "Please specify an output file prefix(-o)" << endl;
			return EXIT_FAILURE;
		}
		
		// The strain taxonomic level is special -- every directory path is assumed to contain
		// a single strain.
		if( ncbi_names_file.empty() ){

			if( ncbi_dir.empty() ){

				cerr << "Please specify the directory to read for ncbi taxonomy files (--ncbi)" << endl;
				return EXIT_FAILURE;
			}

			ncbi_names_file = ncbi_dir + PATH_SEPERATOR + "names.dmp";
		}

		if( ncbi_nodes_file.empty() ){

			if( ncbi_dir.empty() ){

				cerr << "Please specify the directory to read for ncbi taxonomy files (--ncbi)" << endl;
				return EXIT_FAILURE;
			}

			ncbi_nodes_file = ncbi_dir + PATH_SEPERATOR + "nodes.dmp";
		}

		if( ncbi_accession2taxid_file.empty() ){

			if( ncbi_dir.empty() ){

				cerr << "Please specify the directory to read for ncbi taxonomy files (--ncbi)" << endl;
				return EXIT_FAILURE;
			}

			ncbi_accession2taxid_file.push_back(ncbi_dir + PATH_SEPERATOR + "nucl_gb.accession2taxid.gz");
			ncbi_accession2taxid_file.push_back(ncbi_dir + PATH_SEPERATOR + "nucl_wgs.accession2taxid.gz");
		}
		
		// Before we read the NCBI taxonomy files, make sure that we have sequence data to look up
		deque<string> data;
		
		cerr << "Searching for sequence files ... ";
		
		for(int i = optind;i < argc;i++){
			
			// Recursively search argv[i] to find all sequence files
			find_sequence(data, argv[i]);
		}
		
		cerr << "found " << data.size() << " files" << endl;
		
		// Extract the required information from each genBank file (i.e. accession and topology)
		// *before* we parse the taxonomy files (so we don't have to load the entire 
		// accession -> taxid database, which is quite large). Assume that each GenBank file
		// can contain multiple accession
		unordered_multimap<string, SeqInfo> file_to_accession;
		
		cerr << "Parsing sequence file: ";
		
		string update_buffer;
		
		for(deque<string>::const_iterator i = data.begin();i != data.end();++i){
			
			for(string::const_iterator j = update_buffer.begin();j != update_buffer.end();++j){
				cerr << '\b';
			}
			
			for(string::const_iterator j = update_buffer.begin();j != update_buffer.end();++j){
				cerr << ' ';
			}
			
			for(string::const_iterator j = update_buffer.begin();j != update_buffer.end();++j){
				cerr << '\b';
			}
			
			cerr << *i;
			
			update_buffer = *i;
			
			parse_sequence(file_to_accession, *i);
		}
		
		cerr << "\nRead " << file_to_accession.size() << " records in " << data.size() 
			<< " sequence files." << endl;
		
		// Parse the NCBI taxonomy files
		unordered_map<TaxID, string> names;
		deque<TreeNode> tree;
		deque< pair<Accession, TaxID> > accession_to_taxid;
		
		cerr << "Parsing the NCBI taxonomy name file ... ";

		parse_names(names, ncbi_names_file);

		cerr << "done. Found " << names.size() << " {taxid, name} pairs" << endl;
				
		cerr << "Parsing the NCBI taxonomy node file ... ";

		parse_nodes(tree, ncbi_nodes_file);

		cerr << "done." << endl;
				
		// Since there are a large number of accession that we *don't* need, make a
		// list of the accessions that we *do* need
		deque<Accession> target_accessions;
		
		for(unordered_multimap<string, SeqInfo>::const_iterator i = file_to_accession.begin();
			i != file_to_accession.end();++i){
			
			target_accessions.push_back(i->second.accession);
		}
		
		// Sort the target accessions for fast lookup
		sort( target_accessions.begin(), target_accessions.end() );
		
		cerr << "Parsing the NCBI accession to taxid file(s): ";

		for(deque<string>::const_iterator i = ncbi_accession2taxid_file.begin();
			i != ncbi_accession2taxid_file.end();++i){
			
			parse_accession_to_taxid(accession_to_taxid, *i, target_accessions);
		}

		cerr << "\nFound " << accession_to_taxid.size() 
			<< " out of " << target_accessions.size() 
			<< " accession to taxid mappings" << endl;

		if( !ignore_missing_taxid && (target_accessions.size() != accession_to_taxid.size() ) ){
			
			cerr << "The following accessions did *not* have defined taxonomic ids:" << endl;
			
			for(deque<Accession>::const_iterator i = target_accessions.begin();i != target_accessions.end();++i){
				
				deque< pair<Accession, TaxID> >::const_iterator iter = 
					lower_bound( accession_to_taxid.begin(), accession_to_taxid.end(), 
						*i, find_by_accession() );
				
				if( ( iter == accession_to_taxid.end() ) || (iter->first != *i) ){
					
					cout << accession_to_str(*i, true) << endl;
				}
			}
			
			throw __FILE__ ":main: Did not find all of the requested accession to taxid mappings!";
		}
		
		// Write the results sorted by taxonomic label
		for(deque<Taxonomy>::const_iterator taxa_iter = t.begin();taxa_iter != t.end();++taxa_iter){
		
			string taxonomic_level;
			
			switch(*taxa_iter){
				case SUPERKINGDOM:
					taxonomic_level = "superkingdom";
					break;
				case KINGDOM:
					taxonomic_level = "kingdom";
					break;
				case SUBKINGDOM:
					taxonomic_level = "subkingdom";
					break;
				case SUPERPHYLUM:
					taxonomic_level = "superphylum";
					break;
				case PHYLUM:
					taxonomic_level = "phylum";
					break;
				case SUBPHYLUM:
					taxonomic_level = "subphylum";
					break;
				case SUPERCLASS:
					taxonomic_level = "superclass";
					break;
				case CLASS:
					taxonomic_level = "class";
					break;
				case SUBCLASS:
					taxonomic_level = "subclass";
					break;
				case INFRACLASS:
					taxonomic_level = "infraclass";
					break;
				case SUPERORDER:
					taxonomic_level = "superorder";
					break;
				case ORDER:
					taxonomic_level = "order";
					break;
				case PARVORDER:
					taxonomic_level = "parvorder";
					break;
				case SUBORDER:
					taxonomic_level = "suborder";
					break;
				case INFRAORDER:
					taxonomic_level = "infraorder";
					break;
				case SUPERFAMILY:
					taxonomic_level = "superfamily";
					break;
				case FAMILY:
					taxonomic_level = "family";
					break;
				case SUBFAMILY:
					taxonomic_level = "subfamily";
					break;
				case GENUS:
					taxonomic_level = "genus";
					break;
				case SUBGENUS:
					taxonomic_level = "subgenus";
					break;
				case TRIBE:
					taxonomic_level = "tribe";
					break;
				case SUBTRIBE:
					taxonomic_level = "subtribe";
					break;
				case SPECIES:
					taxonomic_level = "species";
					break;
				case SPECIES_GROUP:
					taxonomic_level = "species_group";
					break;
				case SPECIES_SUBGROUP:
					taxonomic_level = "species_subgroup";
					break;
				case SUBSPECIES:
					taxonomic_level = "subspecies";
					break;
				case STRAIN:
					taxonomic_level = "strain";
					break;
				case NO_RANK:
					taxonomic_level = "no_rank";
					break;
				case VARIETAS:
					taxonomic_level = "varietas";
					break;
				case FORMA:
					taxonomic_level = "forma";
					break;
				case UNKNOWN:
					taxonomic_level = "unknown";
					break;
				default:
					throw __FILE__ ":main: Unknown taxonomic level!";
			};
			
			cerr << "Writing gottcha taxa inventory file for: " << taxonomic_level << endl;
			
			gzFile fout = NULL;
			
			if(compress){
				
				const string output_file = output_file_prefix + "." + taxonomic_level + ".txt.gz";
				
				fout = gzopen(output_file.c_str(), "w6");
			}
			else{
				const string output_file = output_file_prefix + "." + taxonomic_level + ".txt";
				
				fout = gzopen(output_file.c_str(), "wT");
			}
			
			if(fout == NULL){

				cerr << "Unable to open output file for writing" << endl;
				return EXIT_FAILURE;
			}

			deque<string> output;

			for(deque<string>::const_iterator i = data.begin();i != data.end();++i){

				// Look up the sequence information associated with this data file and
				// make sure that all of the topology and taxonomy information is the same
				typedef unordered_multimap<string, SeqInfo>::const_iterator I;

				pair<I, I> range = file_to_accession.equal_range(*i);
				
				size_t count = 0;
				
				for(I j = range.first;j != range.second;++j){
					++count;
				}
				
				cerr << "\t" << *i << "Mapping " << count << " sequences" << endl;
				
				string taxonomy_name;

				deque<string> local_output;

				TaxID last_taxid = 0;
				Topology last_topology = UNKNOWN_TOPOLOGY;
				bool single_taxid_and_topology = true;

				for(I j = range.first;j != range.second;++j){

					deque< pair<Accession, TaxID> >::const_iterator iter = 
						lower_bound( accession_to_taxid.begin(), accession_to_taxid.end(),
						j->second.accession, find_by_accession() );

					if( ( iter == accession_to_taxid.end() ) || (iter->first != j->second.accession) ){

						if(ignore_missing_taxid){

							// Skip this record
							continue;
						}

						cerr << "Error in: " << *i << endl;
						cerr << "Accession " << accession_to_str(j->second.accession, true) 
							<< " does not have an associated taxid!" << endl;

						throw __FILE__ ":main: Unable to find taxa id for this accession";
					}

					if(*taxa_iter == STRAIN){

						// Handle the strain as a special case. The taxonomic label is the
						// name of the taxid *for this gi*. However, this approach is far from fool-proof,
						// since while there are examples of multiple genomes of the same strain with the 
						// same taxid, there are also examples of multiple genomes from different strains
						// with the same taxid...
						unordered_map<TaxID, string>::const_iterator name_iter = names.find(iter->second);

						if( name_iter == names.end() ){

							cerr << "** Warning ** Unable to find a species level taxa name for " 
								<< *i << "; taxid = " << iter->second << endl;
							continue;
						}

						taxonomy_name = remove_whitespace(name_iter->second);
					}
					else{
						// Use the sequence taxa id to find the output level taxa id
						const TaxID parent = get_output_taxa(iter->second, tree, *taxa_iter);
						
						if(parent == UNKNOWN){

							cerr << "** Warning ** Unable to find parent taxa id for " 
								<< *i << "; taxa id = " << iter->second << endl;
							continue;
						}

						// Find the name associated with this parent taxa id -- this is the name
						// that will be associated with the sequence file
						unordered_map<TaxID, string>::const_iterator name_iter = names.find(parent);
						
						if( name_iter == names.end() ){

							cerr << "** Warning ** Unable to find a taxa name for " 
								<< *i << "; taxa id = " << iter->second 
								<< "; parent taxa id = " << parent << endl;
							continue;
						}

						taxonomy_name = remove_whitespace(name_iter->second);
					}

					if(j == range.first){

						last_taxid = iter->second;
						last_topology = j->second.topo;
					}
					else{

						if(last_taxid != iter->second){
							single_taxid_and_topology = false;
						}

						if(last_topology != j->second.topo){
							single_taxid_and_topology = false;
						}
					}

					stringstream ssout;

					ssout << taxonomy_name << '\t' << *i << '\t';

					switch(j->second.topo){
						case LINEAR:
						case UNKNOWN_TOPOLOGY: // Assum unknown == linear for now
							ssout << "linear";
							break;
						case CIRCULAR:
							ssout << "circular";
							break;
						default:
							cerr << "Error with file: " << *i << endl;
							throw __FILE__ ": Unknown topology (1)";
					};
					
					ssout << '\t' << accession_to_str(j->second.accession, true) << '\t' << iter->second << '\n';
					
					local_output.push_back( ssout.str() );
				}

				if(single_taxid_and_topology){

					stringstream ssout;

					ssout << taxonomy_name << '\t' << *i << '\t';

					switch(last_topology){
						case LINEAR:
						case UNKNOWN_TOPOLOGY: // Assum unknown == linear for now
							ssout << "linear";
							break;
						case CIRCULAR:
							ssout << "circular";
							break;
						default:
							cerr << "Error with file: " << *i << endl;
							throw __FILE__ ": Unknown topology (1)";
					};

					// A wild card symbol means that *all* accessions in a sequence file map to the same taxonomy and
					// have the same topology
					ssout << '\t' << '*' << '\t' << last_taxid << '\n';

					output.push_back( ssout.str() );
				}
				else{
					output.insert( output.end(), local_output.begin(), local_output.end() );
				}
			}

			sort( output.begin(), output.end() );

			cerr << "Writing " << output.size() << " lines" << endl;

			for(deque<string>::const_iterator i = output.begin();i != output.end();++i){
				
				//fout << *i << endl;
				gzputs( fout, i->c_str() );
			}
			
			gzclose(fout);
		}
		
	}
	catch(const char *error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(const string error){
		cerr << "Caught the error " << error << endl;
		return EXIT_FAILURE;
	}
	catch(...){
		cerr << "Caught an unhandled error" << endl;
		return EXIT_FAILURE;
	}
	
	
	return EXIT_SUCCESS;
}

void parse_accession_to_taxid(deque< pair<Accession, TaxID> > &m_data, const string &m_file, 
	const deque<Accession> &m_target)
{
	gzFile fin = gzopen(m_file.c_str(), "r");
	
	if(fin == NULL){
		throw __FILE__ ":parse_accession_to_taxid: Unable to open accession to taxid file";
	}
	
	// According to the zlib documentation, we can speed up reading by allocating a larger internal
	// buffer for zlib. The gzbuffer function is only availble in more recent versions of zlib
	//if(gzbuffer(fin, 64*1024) != 0){
	//	throw __FILE__ ":parse_accession_to_taxid: Unable to modify zlib buffer";
	//}
	
	const char *status = "/-\\|";
	const size_t status_len = strlen(status);
	const size_t update_every = 1000000;
	size_t status_index = 0;
	
	const int buffer_len = 1024;
	char buffer[buffer_len];
	
	size_t line_number = 0;
	
	try{
		// Skip the first line
		if( !gzgets(fin, buffer, buffer_len) ){
			throw __FILE__ ":parse_accession_to_taxid: Unable to read header";
		}

		// Read data in the form: [accession]\t[accession.version]\t[taxid]
		while( gzgets(fin, buffer, buffer_len) ){

			++line_number;

			if(line_number % update_every == 0){

				cerr << '\b' << status[status_index];

				status_index = (status_index + 1)%status_len;
			}

			TaxID tid = 0;

			char* ptr = buffer;

			while( (*ptr != '\t') && (*ptr != '\0') ){
				++ptr;
			}

			if(*ptr != '\t'){
				throw __FILE__ ":parse_accession_to_taxid: Error reading accession";
			}

			// Skip the tab
			++ptr;

			if(*ptr == '\0'){
				throw __FILE__ ":parse_accession_to_taxid: Error reading delimeter";
			}

			char* accession = ptr;
			size_t accession_len = 0;

			// Read the accession
			while( (*ptr != '\0') && !isspace(*ptr) ){

				++accession_len;
				++ptr;
			}

			if(*ptr != '\t'){
				
				//throw __FILE__ ":parse_accession_to_taxid: Error reading accession.version";
				
				// There is a problem with at least one of the accession to taxid files provided
				// by NCBI (dead_nucl.accession2taxid.gz). Skip this record for now.
				continue;
			}

			const Accession a = str_to_accession(accession, accession + accession_len);

			// Is this accession in the list of accession targets?
			deque<Accession>::const_iterator iter = lower_bound( m_target.begin(),
				m_target.end(), a);

			if( ( iter == m_target.end() ) || (a != *iter) ){
				continue;
			}

			// Skip the tab
			++ptr;

			// Read the taxid
			while( isdigit(*ptr) ){

				tid = tid*10 + (*ptr - '0');
				++ptr;
			}

			// Store this mapping: accession -> tid
			// We can't use an unordered_map due to high memory usage and
			// slow insertion time.
			m_data.push_back( make_pair(a, tid) );
		}
	}
	catch(const char* error){
		
		cerr << "Error parsing " << m_file << " @ line " << line_number << endl;
		cerr << "buffer: " << buffer << endl;
		cerr << "Caught the error: " << error << endl;
		throw error;
	}
	catch(...){
		cerr << "Error parsing " << m_file << " @ line " << line_number << endl;
		cerr << "buffer: " << buffer << endl;
		cerr << "Caught an unhandled error" << endl;
		throw __FILE__ ":parse_accession_to_taxid: Error parsing file";
	}
	
	gzclose(fin);
	
	// Sort the data for fast lookup
	sort( m_data.begin(), m_data.end() );
}

void parse_nodes(deque<TreeNode> &m_tree, const string &m_file)
{
	gzFile fin = gzopen(m_file.c_str(), "r");
	
	if(fin == NULL){
		
		// Try reading a compressed version of this file
		const string compressed = m_file + ".gz";
		
		fin = gzopen(compressed.c_str(), "r");
		
		if(fin == NULL){
			throw __FILE__ ":parse_nodes: Unable to open NCBI node file";
		}
	}
	
	const char delim = '|';
	
	const int buffer_len = 1024;
	char buffer[buffer_len];
	
	unordered_map<string, Taxonomy> ncbi_taxonomy;
	
	ncbi_taxonomy["superkingdom"] = SUPERKINGDOM;
	ncbi_taxonomy["kingdom"] = KINGDOM;
	ncbi_taxonomy["subkingdom"] = SUBKINGDOM;
	ncbi_taxonomy["superphylum"] = SUPERPHYLUM;
	ncbi_taxonomy["phylum"] = PHYLUM;
	ncbi_taxonomy["subphylum"] = SUBPHYLUM;
	ncbi_taxonomy["superclass"] = SUPERCLASS;
	ncbi_taxonomy["class"] = CLASS;
	ncbi_taxonomy["infraclass"] = INFRACLASS;
	ncbi_taxonomy["subclass"] = SUBCLASS;
	ncbi_taxonomy["superorder"] = SUPERORDER;
	ncbi_taxonomy["order"] = ORDER;
	ncbi_taxonomy["parvorder"] = PARVORDER;
	ncbi_taxonomy["suborder"] = SUBORDER;
	ncbi_taxonomy["infraorder"] = INFRAORDER;
	ncbi_taxonomy["superfamily"] = SUPERFAMILY;
	ncbi_taxonomy["subfamily"] = SUBFAMILY;
	ncbi_taxonomy["family"] = FAMILY;
	ncbi_taxonomy["genus"] = GENUS;
	ncbi_taxonomy["subgenus"] = SUBGENUS;
	ncbi_taxonomy["tribe"] = TRIBE;
	ncbi_taxonomy["subtribe"] = SUBTRIBE;
	ncbi_taxonomy["species"] = SPECIES;
	ncbi_taxonomy["species group"] = SPECIES_GROUP;
	ncbi_taxonomy["species subgroup"] = SPECIES_SUBGROUP;
	ncbi_taxonomy["subspecies"] = SUBSPECIES;
	ncbi_taxonomy["strain"] = STRAIN;
	ncbi_taxonomy["no rank"] = NO_RANK;
	ncbi_taxonomy["varietas"] = VARIETAS;
	ncbi_taxonomy["forma"] = FORMA;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		// Split the line according to the format of the nodes.dmp file: 
		//	child taxid | parent taxid | taxa level | ...
		TreeNode node;
		string level;
		string tmp;
		
		char* p = buffer;
		
		// Parse the child taxID
		while( (*p != '\0') && (*p != delim) ){
			
			if( isdigit(*p) ){
				node.child = 10*node.child + (*p - '0');
			}
			
			++p;
		}
		
		// Parse the parent taxID
		if(*p == '\0'){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_names: Unexpected end of line (1)";
		}
		
		++p;
		
		tmp.clear();
		
		while( (*p != '\0') && (*p != delim) ){
			
			if( isdigit(*p) ){
				node.parent = 10*node.parent + (*p - '0');
			}
			
			++p;
		}
		
		// Parse the level
		if(*p == '\0'){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_nodes: Unexpected end of line (2)";
		}
		
		++p;
		
		while( (*p != '\0') && (*p != delim) ){
			
			level.push_back(*p);
			
			++p;
		}
		
		// Remove the flanking white space
		level = tolower( trim_white_space(level) );
		
		// Convert the level to an enumerated taxonomic level
		unordered_map<string, Taxonomy>::const_iterator iter = ncbi_taxonomy.find(level);
		
		if( iter == ncbi_taxonomy.end() ){
			cerr << "Found an undefined taxonomic level: " << level << endl;
		}
		else{
			node.rank = iter->second;
			
			m_tree.push_back(node);
		}
	}
	
	gzclose(fin);
	
	// Sort the nodes in the tree for fast lookup
	sort( m_tree.begin(), m_tree.end() );
}

void parse_names(unordered_map<TaxID, string> &m_names, const string &m_file)
{
	gzFile fin = gzopen(m_file.c_str(), "r");
	
	if(fin == NULL){
		
		// Try reading a compressed version of this file
		const string compressed = m_file + ".gz";
		
		fin = gzopen(compressed.c_str(), "r");
		
		if(fin == NULL){
			throw __FILE__ ":parse_names: Unable to open NCBI name file";
		}
	}
	
	const char delim = '|';
	
	const int buffer_len = 1024;
	char buffer[buffer_len];
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		// Split the line according to the format of the names.dmp file: 
		//	taxid | taxa string | unique taxa string | name class
		TaxID t;
		string taxa;
		string unique_taxa;
		string name_class;
		string tmp;
		
		char* p = buffer;
		
		// Parse the taxID
		while( (*p != '\0') && (*p != delim) ){
			
			if( isdigit(*p) ){
				tmp.push_back(*p);
			}
			
			++p;
		}
		
		t = string_to_int(tmp);
		
		if(t <= 0){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_names: Unable to parse taxID";
		}
		
		// Parse the taxa
		if(*p == '\0'){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_names: Unexpected end of line (1)";
		}
		
		++p;
		
		while( (*p != '\0') && (*p != delim) ){
			
			taxa.push_back(*p);
			
			++p;
		}
		
		// Remove the flanking white space
		taxa = trim_white_space(taxa);
		
		// Parse the unqiue taxa
		if(*p == '\0'){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_names: Unexpected end of line (2)";
		}
		
		++p;
		
		while( (*p != '\0') && (*p != delim) ){
			
			unique_taxa.push_back(*p);
			
			++p;
		}
		
		// Remove the flanking white space
		unique_taxa = trim_white_space(unique_taxa);
		
		// Parse the name class
		if(*p == '\0'){
			
			cerr << "Error with line: " << buffer << endl;
			throw __FILE__ ":parse_names: Unexpected end of line (3)";
		}
		
		++p;
		
		while( (*p != '\0') && (*p != delim) ){
			
			name_class.push_back(*p);
			
			++p;
		}
		
		// Remove the flanking white space
		name_class = trim_white_space(name_class);
		
		// Only use the scientific names
		if(tolower(name_class) == "scientific name"){
			
			// Make sure that this taxid is not a duplicate
			if(m_names.insert( make_pair(t, taxa) ).second == false){
				
				cerr << "Error with line: " << buffer << endl;
				throw __FILE__ ":parse_names: Found a duplicate taxID!";
			}
		}
	}
	
	gzclose(fin);
}

string tolower(const string &m_str)
{
	string ret(m_str);
	
	for(string::iterator i = ret.begin();i != ret.end();++i){
		*i = tolower(*i);
	}
	
	return ret;
}

string trim_white_space(const string &m_string)
{
	const size_t len = m_string.size();
	size_t first;
	size_t last;
	
	for(first = 0;first < len;++first){
		
		if( !isspace(m_string[first]) ){
			break;
		}
	}
	
	for(last = len - 1;last != 0;--last){
		
		if( !isspace(m_string[last]) ){
			break;
		}
	}
	
	if(first < last){
		return m_string.substr(first, last - first + 1);
	}
	
	return string();
}

int string_to_int(const string &m_buffer)
{
	int ret = 0;
	int p = 1;
	
	for(string::const_reverse_iterator i = m_buffer.rbegin();i != m_buffer.rend();++i){
		
		//if( !isdigit(*i) ){
		//	throw __FILE__ ":string_to_int: Invalid character";
		//}
		
		ret += (*i - '0')*p;
		p *= 10;
	}
	
	return ret;
}

void find_sequence(deque<string> &m_data, const string &m_path)
{
	struct stat file_info;

	stat(m_path.c_str(), &file_info);
	
	// Is m_path a directory or a file?
	
	// S_ISREG is true for regular files, NOT directories
	if( S_ISREG(file_info.st_mode) ){
		
		// This is a file. Does it match the list of allowed file
		// extensions?
		const size_t path_len = m_path.size();
		
		for(const FileInfo* f = file_types; f->ext != NULL;++f){
			
			string::size_type loc = m_path.find(f->ext);
			
			if(loc == string::npos){
				continue;
			}
			
			if( loc == path_len - strlen(f->ext) ){
			
				// This is a match!
				m_data.push_back(m_path);
				break;
			}
		}
	}
	else{
		
		// Directory
		DIR *dp;
		struct dirent *dir;

		dp = opendir( m_path.c_str() );

		if(dp == NULL){
		
			cerr << "dir = " << m_path << endl;
			throw __FILE__ ":find_sequence: Unable to open directory";
		}

		while( (dir = readdir(dp)) ){

			// Skip any removed files or directories
			if(dir->d_ino == 0){
				continue;
			}

			// Skip the special directories "." and ".."
			if(strcmp(dir->d_name, ".") == 0){
				continue;
			}

			if(strcmp(dir->d_name, "..") == 0){
				continue;
			}

			find_sequence(m_data, m_path + PATH_SEPERATOR + dir->d_name);
		}

		closedir(dp);
	}
}

// Determine the file type from the extension
FileType get_file_type(const string &m_file)
{
	const size_t len = m_file.size();
	
	for(const FileInfo* f = file_types; f->ext != NULL;++f){
			
		string::size_type loc = m_file.find(f->ext);

		if(loc == string::npos){
			continue;
		}

		if( loc == len - strlen(f->ext) ){

			// This is a match!
			return f->file_type;
		}
	}
	
	return UNKNOWN_FILETYPE;
}

void parse_sequence(unordered_multimap<string, SeqInfo> &m_info, const string &m_file)
{
	switch( get_file_type(m_file) ){
		case GENBANK:
			parse_sequence_gbk(m_info, m_file);
			break;
		case FASTA:
			parse_sequence_fasta(m_info, m_file);
			break;
		default:
			throw __FILE__ ":parse_sequence: Unknown file type!";
	};
}

void parse_sequence_gbk(unordered_multimap<string, SeqInfo> &m_info, const string &m_file)
{
	gzFile fin = gzopen(m_file.c_str(), "r");
	
	if(fin == NULL){
		
		// Try reading a compressed version of this file
		const string compressed = m_file + ".gz";
		
		fin = gzopen(compressed.c_str(), "r");
		
		if(fin == NULL){
			throw __FILE__ ":parse_sequence_gbk: Unable to open GenBank flat file";
		}
	}
	
	const int buffer_len = 1024;
	char buffer[buffer_len];
	
	// We need to extract the following information from a GenBank flat file: 
	//	Topology
	//	Accession
	SeqInfo local;
	
	bool found_topology = false;
	bool found_accession = false;
	size_t record_count = 0;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if( strlen(buffer) == (buffer_len - 1) ){
			
			cerr << "Error parsing " << m_file << " @ line: " << buffer << endl;
			throw __FILE__ ":parse_sequence_gbk: Buffer overflow -> current line is too long!";
		}
		
		// Look for the line that starts with "LOCUS"
		if(strstr(buffer, "LOCUS") == buffer){
			
			if(strstr(buffer, "linear") != NULL){
			
				local.topo = LINEAR;
				found_topology = true;
			}
			else{
				if(strstr(buffer, "circular") != NULL){
				
					local.topo = CIRCULAR;
					found_topology = true;
				}
				else{
					cerr << "Error parsing " << m_file << " @ line: " << buffer << endl;
					throw __FILE__ ":parse_sequence_gbk: Could not find topology information";
				}
			}
		}
		
		// Look for the line that starts with "VERSION"
		if( strstr(buffer, "VERSION") == buffer){
			
			char* ptr = strstr(buffer, "VERSION");
			
			if(ptr == NULL){
				
				cerr << "Error parsing " << m_file << " @ line: " << buffer << endl;
				throw __FILE__ ":parse_sequence_gbk: Could not find VERSION information";
			}
			
			string tmp;
			
			ptr += strlen("VERSION");
			
			// Skip the space between VERSION and the start of the accession
			while( (*ptr != '\0') && isspace(*ptr) ){
				++ptr;
			}
			
			char *start = ptr;
			
			// Find the extent of the accession
			while( (*ptr != '\0') && !isspace(*ptr) ){
				++ptr;
			}
			
			char *stop = ptr;
			
			local.accession = str_to_accession(start, stop);
			found_accession = true;
		}
		
		if(found_topology && found_accession){
			
			m_info.insert( make_pair(m_file, local) );

			++record_count;
			found_topology = false;
			found_accession = false;
		}
	}
	
	gzclose(fin);
	
	if(record_count == 0){
	
		cerr << "Error parsing " << m_file << endl;
		throw __FILE__ ":parse_sequence_gbk: Did not find any records in this file!";
	}	
}

void parse_sequence_fasta(unordered_multimap<string, SeqInfo> &m_info, const string &m_file)
{
	gzFile fin = gzopen(m_file.c_str(), "r");
	
	if(fin == NULL){
		
		// Try reading a compressed version of this file
		const string compressed = m_file + ".gz";
		
		fin = gzopen(compressed.c_str(), "r");
		
		if(fin == NULL){
			throw __FILE__ ":parse_sequence_fasta: Unable to open fasta file";
		}
	}
	
	const int buffer_len = 1048576;
	
	char* buffer = new char [buffer_len];
	
	if(buffer == NULL){
		throw __FILE__ ":parse_sequence_fasta: Unable to allocate buffer";
	}
	
	// Since fasta files do not have topology information, we only need to extract accession
	// information
	SeqInfo local;
	
	local.topo = UNKNOWN_TOPOLOGY;
	
	size_t record_count = 0;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		// Find the fasta deflines
		if(strstr(buffer, ">") == buffer){
		
			if( strlen(buffer) == (buffer_len - 1) ){

				cerr << "Error parsing " << m_file << " @ line: " << buffer << endl;
				throw __FILE__ ":parse_sequence_fasta: Buffer overflow -> current line is too long!";
			}
			
			local.accession = accession_from_fasta_defline(buffer);
			
			m_info.insert( make_pair(m_file, local) );
			
			++record_count;
		}
		
	}
	
	gzclose(fin);
	delete [] buffer;
	
	if(record_count == 0){
	
		cerr << "Error parsing " << m_file << endl;
		throw __FILE__ ":parse_sequence_fasta: Did not find any records in this file!";
	}
}

string remove_whitespace(const string &m_str)
{
	string ret(m_str);
	
	for(string::iterator i = ret.begin();i != ret.end();++i){
		
		if( isspace(*i) ){
			*i = '_';
		}
	}
	
	return ret;
}

deque<string> split(const string &m_path, const char &m_delim)
{
	deque<string> ret;

	string curr;
	
	for(string::const_iterator i = m_path.begin();i != m_path.end();++i){
		
		if(*i == m_delim){
			
			if( !curr.empty() ){
			
				ret.push_back(curr);
				curr.clear();
			}
		}
		else{
			curr.push_back(*i);
		}
	}
	
	if( !curr.empty() ){
			
		ret.push_back(curr);
		curr.clear();
	}

	return ret;
}

TaxID get_output_taxa(const TaxID &m_id, const deque<TreeNode> &m_tree, const Taxonomy &m_level)
{
	TaxID ret = m_id;
	
	// Walk up the taxonomy tree
	while(true){
		
		deque<TreeNode>::const_iterator iter = lower_bound(m_tree.begin(), m_tree.end(), ret);
		
		if( ( iter == m_tree.end() ) || (iter->child != ret) || (ret == iter->parent) ){
			return UNKNOWN;
		}
		
		if(iter->rank == m_level){
			return ret;
		}
		
		ret = iter->parent;
	}
	
	return UNKNOWN;
}
