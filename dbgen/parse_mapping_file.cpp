#include <iostream>
#include <zlib.h>

#include "gottcha_db.h"
#include "deque_set.h"

using namespace std;

struct LocalTaxaRecord
{
	unsigned int group_id;
	string file_name;
	Topology topology;
	TaxaId taxa_id;
	TaxaLevel taxa_level;
	
	// Sort by filename and then by taxa_level
	inline bool operator<(const LocalTaxaRecord &m_rhs) const
	{
		if(file_name == m_rhs.file_name){
			return (taxa_level < m_rhs.taxa_level);
		}
		
		return (file_name < m_rhs.file_name);
	};
};

void parse_id_to_genome(const string &m_filename, const TaxaLevel &m_level,
	unordered_map<string, unsigned int> &m_group_name_to_id,
	unsigned int &m_next_group_index, deque<LocalTaxaRecord> &m_buffer, 
	const Topology &m_default_topo);
string trim_white_space(const string &m_string);

void parse_mapping_files(vector< pair<SequenceFileInfo, string> > &m_file_info, 
	deque<unsigned int> &m_groups_to_squash,
	const Options &m_opt, UpdateInfo &m_info)
{
	deque<LocalTaxaRecord> local;
	
	unordered_map<string, unsigned int> group_name_to_id; // Used to assign unique group name ids
	unsigned int next_group_index = 1; // Index groups starting from 1 (0 is a special value)
	
	// Parse the mapping file(s) for each taxa level
	for(unordered_map<TaxaLevel, string>::const_iterator i = m_opt.taxa_level_mapping_file.begin();
		i != m_opt.taxa_level_mapping_file.end();++i){

		if(m_opt.verbose >= Options::INFORMATIVE){

			m_info << "Reading " << taxa_level_name(i->first) 
				<< "-level file mapping information";
				
			m_info.flush();
		}
		
		parse_id_to_genome(i->second, i->first, group_name_to_id, 
			next_group_index, local,
			m_opt.default_topology);
	}
	
	// Save the group ids to squash
	for(deque<string>::const_iterator i = m_opt.groups_to_squash.begin();
		i != m_opt.groups_to_squash.end();++i){
		
		unordered_map<string, unsigned int>::const_iterator iter = group_name_to_id.find(*i);
		
		if( iter == group_name_to_id.end() ){
			
			cerr << "Did not find the group name to squash, " << *i << ", in the list of actual groups" << endl;			
			throw __FILE__ ":parse_mapping_files: Unable to find requested group to squash";
		}
		
		m_groups_to_squash.push_back(iter->second);
	}
	
	make_set(m_groups_to_squash);
	
	// Sort the filename mapping records so we can get a count of the number of unique filenames
	SORT( local.begin(), local.end() );
	
	// We will only index the unique filenames
	string last_filename;
	
	size_t num_unique_files = 0;
	
	for(deque<LocalTaxaRecord>::const_iterator i = local.begin();i != local.end();++i){
		
		if( i->file_name != last_filename ){
		
			++num_unique_files;
			last_filename = i->file_name;
		}
	}
		
	m_file_info.resize(num_unique_files);
	
	last_filename = "";
	size_t file_name_index = 0;
	
	for(deque<LocalTaxaRecord>::const_iterator i = local.begin();i != local.end();++i){
		
		if( i == local.begin() ){
			last_filename = i->file_name;
		}
		else{
			if( i->file_name != last_filename ){

				++file_name_index;
				last_filename = i->file_name;
			}
		}
		
		m_file_info[file_name_index].second = i->file_name;
		
		SequenceFileInfo &info = m_file_info[file_name_index].first;
		
		info.group_id[i->taxa_level] = i->group_id;
		
		if(info.topo == UNKNOWN_TOPOLOGY){
		
			info.topo = i->topology;
			info.taxa_id = i->taxa_id;
		}
		else{
			// Make sure that the topology and taxa_id information is consistent
			// across all taxonomic levels
			if(info.topo != i->topology){
				
				cerr << "\nError @ " << i->file_name << endl;
				cerr << "\tPrevious topology (" << int(info.topo) 
					<< "); current topology ("
					<< int(i->topology) << ")" << endl;
				throw __FILE__ ":parse_mapping_files: Inconsistent file topology";
			}
			
			if(info.taxa_id != i->taxa_id){
				
				cerr << "\nError @ " << i->file_name << endl;
				cerr << "\tPrevious taxa id (" << info.taxa_id
					<< "); current taxa id ("
					<< i->taxa_id << ")" << endl;
				throw __FILE__ ":parse_mapping_files: Inconsistent file taxa id";
			}
		}
	}
}

void parse_id_to_genome(const string &m_filename, const TaxaLevel &m_level,
	unordered_map<string, unsigned int> &m_group_name_to_id,
	unsigned int &m_next_group_index, deque<LocalTaxaRecord> &m_buffer, 
	const Topology &m_default_topo)
{
	// Allow compressed id to genome files
	gzFile fin = gzopen( m_filename.c_str(), "r");
	
	if(fin == NULL){
	
		cerr << "Unable to open taxa id to genome file (" << m_filename << ") for reading" << endl;
		throw __FILE__ ":parse_id_to_genome: Unable to open file";
	}
	
	// The group id, path and topology fields must be non-space delimited, since a 
	// space is a valid component of a filename. Use tabs for now...
	const char delim = '\t';

	size_t line_number = 0;
		
	const size_t buffer_len = 65536;
	char buffer[buffer_len];
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		try{
			const string line(buffer);
			
			string group_name;
			string file_name;
			string topology;
			string taxa_id;
			
			++line_number;

			if( line.size() == (buffer_len - 1) ){
				throw __FILE__ ":parse_id_to_genome: Buffer to small to read entire line!";
			}
		
			string::const_iterator i = line.begin();

			// Skip any leading white space
			for(;( i != line.end() ) && isspace(*i);++i){
			}

			// Skip blank lines and lines that start with a comment
			if( ( i == line.end() ) || (*i == '#') ){
				continue;
			}
			
			// Read the group id -- break on the delimeter, *not* white space!
			for(;( i != line.end() ) && (*i != delim);++i){
				group_name.push_back(*i);
			}

			// Remove any flanking white space from the group id
			group_name = trim_white_space(group_name);
			
			// Skip the separating white space
			for(;( i != line.end() ) && isspace(*i);++i){
			}

			// Read the file path -- break on the delimeter, *not* white space!
			for(;( i != line.end() ) && (*i != delim);++i){
				file_name.push_back(*i);
			}

			// Skip the separating white space
			for(;( i != line.end() ) && isspace(*i);++i){
			}

			// Read the topology
			for(;( i != line.end() ) && !isspace(*i);++i){
				topology.push_back( toupper(*i) );
			}

			// Allow the user to specify a default topology that overrides
			// the topology specified in the input file
			switch(m_default_topo){
				case CIRCULAR:
					topology = "CIRCULAR";
					break;
				case LINEAR:
					topology = "LINEAR";
					break;
				case UNKNOWN_TOPOLOGY:
					// Do nothing (the input file specifies the topology in this case)
					break;
			};
			
			// Skip the separating white space
			for(;( i != line.end() ) && isspace(*i);++i){
			}

			// Read the taxa id
			for(;( i != line.end() ) && !isspace(*i);++i){
				taxa_id.push_back(*i);
			}
			
			if( group_name.empty() || file_name.empty() || 
				topology.empty() || taxa_id.empty() ){

				cerr << "Error parsing taxa id to genome file at line number " << line_number << endl;

				if( group_name.empty() ){
					cerr << "\tMissing group taxa id" << endl;
				}

				if( file_name.empty() ){
					cerr << "\tMissing path to genome file" << endl;
				}

				if( topology.empty() ){
					cerr << "\tMissing genome topology (remember to use a tab, not spaces)" << endl;
				}
				
				if( taxa_id.empty() ){
					cerr << "\tMissing taxa id" << endl;
				}
				
				throw __FILE__ ":parse_id_to_genome: I/O error";
			}

			Topology local_topo = UNKNOWN_TOPOLOGY;
			
			if(topology == "LINEAR"){
				local_topo = LINEAR;
			}
			else{
				if(topology == "CIRCULAR"){
					local_topo = CIRCULAR;
				}
				else{
					throw __FILE__ ":parse_id_to_genome: Unknown topology";
				}
			}
			
			// Is this a new group?
			unordered_map<string, unsigned int>::const_iterator group_iter =
				m_group_name_to_id.find(group_name);
			
			unsigned int group_id = 0xFFFFFFFF; // An invalid group id
			
			if( group_iter == m_group_name_to_id.end() ){

				// This is a new group
				group_id = m_next_group_index;

				m_group_name_to_id.insert( make_pair(group_name, group_id) );

				++m_next_group_index;
			}
			else{
				// This is an existing group
				group_id = group_iter->second;
			}
			
			// Store this record
			m_buffer.push_back( LocalTaxaRecord() );
			
			LocalTaxaRecord &ref = m_buffer.back();
			
			ref.group_id = group_id;
			ref.file_name = file_name;
			ref.topology = local_topo;
			ref.taxa_id = taxa_id;
			ref.taxa_level = m_level;
		}
		catch(const char* error){
			
			cerr << "Error parsing line #" << line_number << " in file: " << m_filename << endl;
			gzclose(fin);
			throw error;
		}
		catch(...){
			
			cerr << "Error parsing line #" << line_number << " in file: " << m_filename << endl;
			gzclose(fin);
			throw __FILE__ ":Options::parse_id_to_genome: Caught an unhandled error!";
		}
	}

	gzclose(fin);
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
	
	if(first <= last){
		return m_string.substr(first, last - first + 1);
	}
	
	return string();
}
