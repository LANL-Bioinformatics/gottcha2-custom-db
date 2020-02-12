#include <fstream>
#include <algorithm>
#include <list>
#include <unordered_set>
#include <iostream>

#include <mpi.h>
#include <string.h>

// For testing and creating directories
#include <sys/types.h>
#include <sys/stat.h>
#include <dirent.h>
#include <ctype.h>

#include "gottcha_db.h"
#include "deque_set.h"

using namespace std;

struct Capacity
{

	int rank;
	size_t size;
	
	Capacity()
	{
		rank = 0;
		size = 0;
	};
	
	Capacity(const int &m_rank, const size_t &m_size) :
		rank(m_rank), size(m_size)
	{

	};
	
	inline bool operator<(const Capacity &m_rhs) const
	{
		// Change the default ordering when sorting (such that
		// the smallest element is sorted into the first position).
		return size > m_rhs.size;
	};
};

bool extension_match(const string &m_path, const string &m_ext);
string compute_lcp(const vector< pair<SequenceFileInfo, string> > &m_file_info);
pair<unsigned char*, size_t> load_sequence_buffer(const string &m_filename,
	const SequenceFileInfo &m_info, const size_t &m_word_size);
void send_sequence_file(const int &m_target_db_rank, const size_t &m_file_index, 
	const string &m_file_name, const SequenceFileInfo &m_file_info);

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

void load_sequence(deque< pair<string, SequenceZ*> > &m_sequence, 
	const string &m_filename, const SequenceFileInfo &m_file_info,
	const size_t& m_word_size)
{
	switch( GetFileType(m_filename) ){
		case FASTA:
		
			parse_fasta(m_filename, m_sequence, m_file_info, m_word_size);
			break;
		case GENBANK:
			
			parse_genbank(m_filename, m_sequence, m_file_info, m_word_size);
			break;
		default:
			cerr << "m_filename = " << m_filename << endl;
			throw __FILE__ ":parse_sequence: Unknown file type!";
	};
}
	
bool extension_match(const string &m_path, const string &m_ext)
{
	const size_t ext_len = m_ext.size();
	
	if( m_path.size() < ext_len ){
		return false;
	}
	
	string::const_reverse_iterator path_iter = m_path.rbegin();
	
	for(string::const_reverse_iterator ext_iter = m_ext.rbegin();
		ext_iter != m_ext.rend();++ext_iter, ++path_iter){
		
		if(*ext_iter != *path_iter){
			return false;
		}
	}
	
	return true;
}

string remove_extension(const string &m_name, const string &m_ext)
{
	string::size_type loc = m_name.rfind(m_ext);
	
	if( (loc != string::npos) && ( loc == ( m_name.size() - m_ext.size() ) ) ){
		
		// *Remove* the extension
		return m_name.substr(0, loc);
	}
	
	return m_name;
}

FileType GetFileType(const string &m_path)
{
	// Remove a file compression extension (if any)
	const string path = remove_extension(m_path, ".gz");
	
	// For now, use simple file extension matching
	if( extension_match(path, ".fna") || 
	    extension_match(path, ".fa") ||
	    extension_match(path, ".fasta") ){
		return FASTA;
	}
		
	if( extension_match(path, ".gbk") ||
	    extension_match(path, ".gb") ||
	    extension_match(path, ".gbff") ){
		return GENBANK;
	}
	
	return UNKNOWN_FILETYPE;
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
			
			// Save the delimeter too!
			curr.push_back(*i);
			ret.push_back(curr);
			curr.clear();
		}
		else{
			curr.push_back(*i);
		}
	}
	
	if( !curr.empty() ){
		ret.push_back(curr);
	}
	
	return ret;
}

// Return the longest common prefix of the input files (if a new root directory has
// been specified. Otherwise return an empty string.
string load_sequence_database(
	unordered_map<size_t, pair<unsigned char*, size_t> > &m_db,
	vector<int> &m_sequence_location,
	deque<unsigned int> &m_groups_to_squash,
	UpdateInfo &m_progress, const Options &m_opt)
{
	string longest_common_prefix;
	
	const double profile_load = MPI_Wtime();
	
	vector< pair<SequenceFileInfo, string> > file_info;
	size_t num_seq_files = 0;
	
	// Only rank 0 needs to read the mapping files
	if(mpi_rank == 0){
	
		parse_mapping_files(file_info, m_groups_to_squash, 
			m_opt, m_progress);
		
		// Sort the input files according to taxonomic groups. Grouping multiple sequences
		// that share taxonomic groups is expected to improve performance
		SORT( file_info.begin(), file_info.end() );
		
		num_seq_files = file_info.size();
		
		// If the user has specified a new root directory, identify
		// the longest common prefix for all of the input files (which will
		// be stripped off and replaced with the new root directory when the
		// output files are written).
		if( !m_opt.output_root_dir.empty() ){
		
			longest_common_prefix = compute_lcp(file_info);
		
			// Since we will be writing to a new root directory, we will potentially
			// need to create a new directory structure to store the files. Have the
			// rank 0 process perform this task to avoid any potential race conditions
			// between ranks if they had to create output directories.
			unordered_set<string> dir_tree;
			
			for(vector< pair<SequenceFileInfo, string> >::const_iterator i = file_info.begin();
				i != file_info.end();++i){
				
				// Skip files were every taxonomic level has been squashed
				bool all_squashed = true;
				
				for(unsigned int level = 0;(level < NUM_TAXA_LEVEL) && all_squashed;++level){
					
					// Keep testing for squashed groups until we find a group that is *not*
					// squashed
					all_squashed = set_contains(m_groups_to_squash, i->first.group_id[level]);
				}
				
				if(all_squashed){
					
					// Since all taxonomic levels for this file have been squashed, we do
					// not need to make an output file for it.
					continue;
				}
				
				const string filename = format_output_filename(i->second, 
					"dummy", m_opt.output_root_dir, 
					longest_common_prefix, m_opt.compress_output);
				
				const size_t len = filename.size();
	
				// Since runs of multiple PATH_SEPARATORs are legal, we only 
				// want to call create_directory when there is an
				// intervening character between adjacent PATH_SEPARATOR
				// symbols.
				size_t last_separator = 0;

				for(size_t j = 0;j < len;++j){

					if(filename[j] == PATH_SEPARATOR){

						if( (j - last_separator) > 1 ){

							// Create a directory if it does not already exist
							dir_tree.insert( filename.substr(0, j) );							
						}

						last_separator = j;
					}
				}
			}
			
			vector<string> dir_to_create( dir_tree.begin(), dir_tree.end() );
			
			dir_tree.clear();
			
			// Sort the directories in ascending order to make sure we create the parent
			// (i.e. shorter) directories first.
			SORT( dir_to_create.begin(), dir_to_create.end() );
			
			const size_t num_dir = dir_to_create.size();
			const size_t update_every = max(size_t(1), num_dir/100);
			size_t dir_count = 0;
			
			for(vector<string>::const_iterator i = dir_to_create.begin();
				i != dir_to_create.end();++i, ++dir_count){
				
				create_directory(*i);
				
				// Keep the user informed
				if( (dir_count%update_every == 0) && 
				    (m_opt.verbose >= Options::INFORMATIVE) ){

					m_progress << "Building new directory tree ("
						<< (100.0*dir_count)/num_dir << "%): "
						<< truncate_filename(*i, 75);
					m_progress.flush();
				}

			}
			
			// Keep the user informed
			if(m_opt.verbose >= Options::INFORMATIVE){

				m_progress << dir_to_create.size() << " directories created";
				m_progress.flush();
			}
		}
		
		// Keep the user informed
		if(m_opt.verbose >= Options::INFORMATIVE){
		
			m_progress << "Found " << file_info.size() << " sequence files to process" << endl;
			
			if( !m_opt.output_root_dir.empty() ){
				m_progress << "Output files will replace prefix \"" << longest_common_prefix 
					<< "\" with \"" << m_opt.output_root_dir << '"' << endl;
			}
		}		
	}
	
	// Share the number of files with the other compute ranks
	MPI_Bcast(&num_seq_files, sizeof(num_seq_files), MPI_BYTE, 0, MPI_COMM_WORLD);
	
	// Share the m_groups_to_squash with the other ranks
	broadcast(m_groups_to_squash, mpi_rank, 0);
	
	// Share the longest common prefix of all input filenames
	broadcast(longest_common_prefix, mpi_rank, 0);
	
	// Make room to store the location of each sequence file
	m_sequence_location.resize(num_seq_files);

	// Compute rank 0 is in charge of orchestrating the distribution of sequence files amoung
	// the different database ranks	
	if(mpi_rank == 0){
		
		// All of the ranks need to send their hostname to rank 0. To avoid overtaxing
		// the disk I/O subsystem, we will only allocate sequence data to at most N ranks
		// per physical host
		const size_t max_ranks_per_host = 2; // Should be a user-defined option!
		unordered_map<string, size_t> ranks_per_host;
		
		// Maintain a priority queue of database ranks, such that the rank with the smallest amount of sequence data
		// is preferentially selected for the next file
		priority_queue<Capacity> db_ready;
		
		// Skip adding rank 0 for now (it will be added later, as a special case)
		for(int i = 1;i < mpi_numtasks;++i){
			
			char hostname[MPI_MAX_PROCESSOR_NAME];
			
			if(MPI_Recv(hostname, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, i, 
				GOTTCHA_HOSTNAME, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){
		
				throw __FILE__ ":load_sequence_database: Error receiving rank hostname";
			}
			
			size_t &rank_count = ranks_per_host[hostname];
					
			++rank_count;

			if(rank_count <= max_ranks_per_host){
				db_ready.push( Capacity(i, 0) );
			}
		}

		vector<size_t> db_size(mpi_numtasks);
		
		MPI_Status status;
		
		// Update the user every 1% of the number of target files
		const size_t update_every = max( size_t(1), size_t(0.01*num_seq_files) );
		
		// Keep track of the number of ranks that are currenting working on loading data
		// (there are some large genomes that make take a lot of time).
		int num_working = 0;
		
        	for(size_t i = 0;i < num_seq_files;++i){

			int target_rank = -1;

			while(true){

				// Are there any ranks that are finished, and ready to receive additional
				// sequence files?
				int db_ret = true;

				while(db_ret){

					if(MPI_Iprobe(MPI_ANY_SOURCE, GOTTCHA_FILE_SIZE, MPI_COMM_WORLD,
						&db_ret, &status) != MPI_SUCCESS){

						throw __FILE__ ":load_sequence_database: Error in MPI_Iprobe(GOTTCHA_FILE_SIZE)";
					}

					if(db_ret){

						size_t file_size;

						if(MPI_Recv(&file_size, sizeof(size_t), MPI_BYTE, status.MPI_SOURCE, 
							GOTTCHA_FILE_SIZE, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

							throw __FILE__ ":load_sequence_database: Error receiving file size";
						}

						db_size[status.MPI_SOURCE] += file_size;

						db_ready.push( Capacity(status.MPI_SOURCE, db_size[status.MPI_SOURCE]) );
						
						// This rank is no longer working
						--num_working;
					}
				}

				// Rather than waiting for rank 0 to load a file, only allow rank 0 to load
				// when it is the only rank.
				if(mpi_numtasks == 1){
                		    db_ready.push( Capacity(0, db_size[0]) );;
                		}
                
                		if( db_ready.empty() ){

                        		// Wait for a rank to become free
                        		continue;
				}

				target_rank = db_ready.top().rank;
				db_ready.pop();

				break;
			}
            
			// Send the information about the current file to the target_db_rank
			if(target_rank == 0){
			
				m_db[i] = load_sequence_buffer(file_info[i].second, file_info[i].first, 
					m_opt.word_size);
				
				db_size[0] += m_db[i].second;
			}
			else{
				send_sequence_file(target_rank, i, file_info[i].second, file_info[i].first);
				
				// This rank is working
				++num_working;
			}

			if( (m_opt.verbose >= Options::INFORMATIVE) && (i%update_every == 0) ){

				m_progress << "Loading " << truncate_filename(file_info[i].second, 75) << " on rank [" 
					<< target_rank << "]: " << (i*100.0)/num_seq_files << "%";
				m_progress.flush();
			}

			m_sequence_location[i] = target_rank;
		}
		
		// Wait for all of the ranks to finish loading so we can report an accurate count
		// of file sizes
		while(num_working != 0){

			size_t file_size;

			if(MPI_Recv(&file_size, sizeof(size_t), MPI_BYTE, MPI_ANY_SOURCE, 
				GOTTCHA_FILE_SIZE, MPI_COMM_WORLD, &status) != MPI_SUCCESS){

				throw __FILE__ ":load_sequence_database: Error receiving file size";
			}

			db_size[status.MPI_SOURCE] += file_size;

			// This rank is no longer working
			--num_working;
		}

		// Tell the other ranks that we are done with sequence distribution by sending a zero
		// length buffer
		for(int i = 1;i < mpi_numtasks;++i){
		
			size_t zero = 0;
			
			if(MPI_Send(&zero, sizeof(size_t), MPI_BYTE, i, 
				GOTTCHA_FILE_LOAD, MPI_COMM_WORLD) != MPI_SUCCESS){

				throw __FILE__ ":load_sequence_database: Unable to zero length buffer size";
			}
		}
		
		if(m_opt.verbose >= Options::INFORMATIVE){
	
			// Collect the amount of data stored on each database rank
			size_t total_db_size = 0;

			// Rank 0 only loads sequence files when it is the only rank
			for(int i = (mpi_numtasks == 1) ? 0 : 1;i < mpi_numtasks;++i){

				total_db_size += db_size[i];

				if(db_size[i] != 0){
				
					m_progress << "\tRank[" << i << "] loaded " << double(db_size[i])/GB
						<< " GB" << endl;
					m_progress.flush();
				}
			}

			m_progress << "Sequence loading is complete: " << double(total_db_size)/GB 
				<< " GB in " << MPI_Wtime() - profile_load << " sec" << endl;
			m_progress.flush();
		}
	}
	else{ // mpi_rank != 0
		
		// Send our hostname to rank 0
		char hostname[MPI_MAX_PROCESSOR_NAME];
		int hostname_len = 0;
		
		MPI_Get_processor_name(hostname, &hostname_len);
		
		if(MPI_Send(hostname, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, 0, 
			GOTTCHA_HOSTNAME, MPI_COMM_WORLD) != MPI_SUCCESS){
		
			throw __FILE__ ":load_sequence_database: Unable to send hostname";
		}
			
		while(true){
			
			size_t buffer_size = 0;
			
			if(MPI_Recv(&buffer_size, sizeof(size_t), MPI_BYTE, 0, 
				GOTTCHA_FILE_LOAD, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

				throw __FILE__ ":load_sequence_database: Error receiving buffer size";
			}
			
			// Rank 0 will tell us when there are no more sequences to load by sending
			// a zero length buffer
			if(buffer_size == 0){
				break;
			}
			
			unsigned char *buffer = new unsigned char [buffer_size];
			
			if(buffer == NULL){
				throw __FILE__ ":load_sequence_database: Error receiving buffer";
			}
			
			if(MPI_Recv(buffer, buffer_size, MPI_BYTE, 0, 
				GOTTCHA_FILE_LOAD, MPI_COMM_WORLD, MPI_STATUS_IGNORE) != MPI_SUCCESS){

				throw __FILE__ ":load_sequence_database: Error receiving buffer";
			}
			
			size_t file_index;
	
			unsigned char* ptr = buffer;

			ptr = mpi_unpack(ptr, file_index);

			string name;
			SequenceFileInfo info;

			ptr = mpi_unpack(ptr, name);
			ptr = mpi_unpack(ptr, info);
			
			delete [] buffer;
			
			m_db[file_index] = load_sequence_buffer(name, info, m_opt.word_size);
			
			// Let rank 0 know we're done by sending the size of sequence file we just loaded
			if(MPI_Send(&(m_db[file_index].second), sizeof(size_t), MPI_BYTE, 0, 
				GOTTCHA_FILE_SIZE, MPI_COMM_WORLD) != MPI_SUCCESS){

				throw __FILE__ ":load_sequence_database: Unable to send file size";
			}
		}
	}
	
	// Share the location of all sequence files will all ranks
	broadcast(m_sequence_location, mpi_rank, 0);
	
	return longest_common_prefix;
}

// Compute the longerst common prefix for the input files
string compute_lcp(const vector< pair<SequenceFileInfo, string> > &m_file_info)
{
	vector< pair<SequenceFileInfo, string> >::const_iterator i = m_file_info.begin();
	
	if( i == m_file_info.end() ){
		return string();
	}
	
	string ret = i->second;
	++i;
	
	size_t prefix_len = ret.size();
	
	for(;i != m_file_info.end();++i){
		
		const size_t len = min( prefix_len, i->second.size() );
		
		size_t j = 0;
		
		for(; (j < len) && (ret[j] == i->second[j]) ;++j){
		}
		
		if(prefix_len != j){
		
			ret = ret.substr(0, j);
			prefix_len = j;
			
			if(prefix_len == 0){
				break;
			}
		}
	}
	
	// Chop any characters after the last path separator
	const string::size_type loc = ret.find_last_of(PATH_SEPARATOR);
	
	if(loc != string::npos){
		ret = ret.substr(0, loc);
	}
	
	return ret;
}

string format_output_filename(const string &m_input_filename, const string &m_prefix, 
	const string &m_root_dir, const string &m_input_lcp, const bool &m_compress)
{
	if( m_prefix.empty() ){
		throw __FILE__ ":Options::format_output: Empty output prefix!";
	}
	
	// Initialize the output path from the input path
	deque<string> path = split(m_input_filename, PATH_SEPARATOR);

	if( path.empty() ){
		throw __FILE__ ":format_output_filename: Unable to process output path!";
	}

	string &filename = path.back();

	// Remove any existing prefix
	filename = remove_extension(filename, ".gz");

	filename = remove_extension(filename, ".fna");
	filename = remove_extension(filename, ".fa");
	filename = remove_extension(filename, ".gb");
	filename = remove_extension(filename, ".gbk");
	filename = remove_extension(filename, ".gbff");

	if(m_compress){
		filename = m_prefix + filename + ".fna.gz";
	}
	else{
		filename = m_prefix + filename + ".fna";
	}

	// Replace the old output filename with the new (correct formatted one) one
	string final_output;

	if( m_root_dir.empty() ){

		// Note that the path deque *includes* the path delimiter
		for(deque<string>::const_iterator i = path.begin();i != path.end();++i){
			final_output += *i;
		}
	}
	else{		
		// Replace the longest common input prefix with the new root directory
		const deque<string> input_path = split(m_input_lcp, PATH_SEPARATOR);
		
		deque<string>::const_iterator j = input_path.begin();
		deque<string>::const_iterator i = path.begin();
		
		for(;i != path.end();++i, ++j){
			
			if( ( j == input_path.end() ) || (*i != *j) ){
				break;
			}
		}
		
		final_output = m_root_dir;
		
		for(;i != path.end();++i){
			final_output += *i;
		}
	}

	return final_output;
}

void create_complete_path(const string &m_filename)
{
	const size_t len = m_filename.size();
	
	// Since runs of multiple PATH_SEPARATORs are legal, we only 
	// want to call create_directory when there is an
	// intervening character between adjacent PATH_SEPARATOR
	// symbols.
	size_t last_separator = 0;
	
	for(size_t i = 0;i < len;++i){
		
		if(m_filename[i] == PATH_SEPARATOR){
			
			if( (i - last_separator) > 1 ){
			
				// Create a directory if it does not already exist
				create_directory( m_filename.substr(0, i) );
			}
			
			last_separator = i;
		}
	}
}

void create_directory(const string &m_path)
{
	// Only read files (need to exclude directories)
	struct stat path_info;

	if(stat(m_path.c_str(), &path_info) == 0){
		
		// A file or directory with this name already exists
		if( !S_ISDIR(path_info.st_mode) ){
			throw __FILE__ ":create_directory: Cannot replace a file with a directory!";
		}
		
		return;
	}

	// Create this directory
	if(mkdir(m_path.c_str(), 0777) != 0){
		throw __FILE__ ":create_directory: Unable to create directory";
	}
}

string truncate_filename(const string &m_name, const size_t &m_max_len)
{
	const size_t len = m_name.size();
	
	if(len <= m_max_len){
		return m_name;
	}
	
	return string("...") + m_name.substr(len - m_max_len, m_max_len);
	
}

pair<unsigned char*, size_t> load_sequence_buffer(const string &m_filename,
	const SequenceFileInfo &m_info, const size_t &m_word_size)
{
	SequenceFile local;
	
	local.name = m_filename;
	local.info = m_info;
	
	load_sequence(local.seq, local.name, local.info, m_word_size);
	
	const size_t buffer_size = mpi_size(local);
	
	if(buffer_size > INT_MAX){
		throw __FILE__ ":load_sequence_buffer: 32 bit MPI buffer overflow!";
	}
	
	// This memory must be deallocated by the calling function
	unsigned char *buffer = new unsigned char [buffer_size];
	
	if(buffer == NULL){
		throw __FILE__ ":load_sequence_buffer: Unable to allocate buffer";
	}
	
	mpi_pack(buffer, local);
	
	// We need to manually deallocate the sequence data in the local SequenceFile
	local.clear();
	
	return make_pair(buffer, buffer_size);
}

void send_sequence_file(const int &m_target_db_rank, const size_t &m_file_index, 
	const string &m_file_name, const SequenceFileInfo &m_file_info)
{

	const size_t buffer_size = sizeof(m_file_index) + mpi_size(m_file_name) + mpi_size(m_file_info);
	
	if(MPI_Send(&buffer_size, sizeof(size_t), MPI_BYTE, m_target_db_rank, 
		GOTTCHA_FILE_LOAD, MPI_COMM_WORLD) != MPI_SUCCESS){
		
		throw __FILE__ ":send_sequence_file: Unable to send buffer size";
	}
	
	unsigned char* buffer = new unsigned char[buffer_size];
	
	if(buffer == NULL){
		throw __FILE__ ":send_sequence_file: Unable to allocate send buffer";
	}
	
	unsigned char* ptr = buffer;
	
	ptr = mpi_pack(ptr, m_file_index);
	ptr = mpi_pack(ptr, m_file_name);
	ptr = mpi_pack(ptr, m_file_info);
	
	if(MPI_Send(buffer, buffer_size, MPI_BYTE, m_target_db_rank, 
		GOTTCHA_FILE_LOAD, MPI_COMM_WORLD) != MPI_SUCCESS){
		
		throw __FILE__ ":send_sequence_file: Unable to send buffer";
	}
	
	delete [] buffer;
}

void request_file(const size_t &m_file_index, const int &m_src_rank, SequenceFile &m_background, 
	const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db)
{

	unsigned char *buffer = NULL;
	
	// Every rank needs to load the current target which is stored by m_target_rank
	if(mpi_rank == m_src_rank){
	
		unordered_map<size_t, pair<unsigned char*, size_t> >::const_iterator iter = 
			m_db.find(m_file_index);
		
		if( iter == m_db.end() ){
			throw __FILE__ ":request_file: Unable to find requested sequence file";
		}
		
		buffer = iter->second.first;
		
		MPI_Bcast((void*)&(iter->second.second), sizeof(size_t), MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
		MPI_Bcast(buffer, iter->second.second, MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
	}
	else{
		
		size_t buffer_size = 0;
		
		MPI_Bcast(&buffer_size, sizeof(size_t), MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
		
		buffer = new unsigned char [buffer_size];
		
		if(buffer == NULL){
			throw __FILE__ ":request_file: Unable to allocate buffer";
		}
		
		MPI_Bcast(buffer, buffer_size, MPI_BYTE, m_src_rank, MPI_COMM_WORLD);
		
	}
	
	mpi_unpack(buffer, m_background);
		
	if(mpi_rank != m_src_rank){
		delete [] buffer;
	}
}

size_t per_rank_database_size(const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db)
{
	// Compute the size of the database stored on each *host* and then divide this size by the
	// number ranks on the host. Note that the local_database_size can be zero for a given rank.
	size_t local_database_size = m_db.size()*2*sizeof(size_t); 
		
	for(unordered_map<size_t, pair<unsigned char*, size_t> >::const_iterator i = m_db.begin();
		i != m_db.end();++i){

		local_database_size += i->second.second;
	}
	
	char hostname[MPI_MAX_PROCESSOR_NAME];
	int hostname_len = 0;
		
	MPI_Get_processor_name(hostname, &hostname_len);
	
	size_t host_database_size = 0;
	size_t ranks_per_host = 0;
	
	for(int i = 0;i < mpi_numtasks;++i){
		
		char host[MPI_MAX_PROCESSOR_NAME];
		size_t db_size = 0;
		
		strcpy(host, hostname);
		db_size = local_database_size;
		
		MPI_Bcast(host, MPI_MAX_PROCESSOR_NAME, MPI_BYTE, i, MPI_COMM_WORLD);
		MPI_Bcast(&db_size, sizeof(size_t), MPI_BYTE, i, MPI_COMM_WORLD);
		
		if(strcmp(host, hostname) == 0){
		
			host_database_size += db_size;
			++ranks_per_host;
		}
	}
	
	return host_database_size/ranks_per_host;
}

