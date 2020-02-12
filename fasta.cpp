#include <fstream>
#include <vector>
#include <iostream>
#include <algorithm>
#include <string.h>

#include <zlib.h>

#include "gottcha_db.h"
#include "deque_set.h"

using namespace std;

void parse_fasta(const string &m_filename, 
	deque< pair<std::string, SequenceZ*> > &m_sequence, 
	const SequenceFileInfo &m_file_info, const size_t& m_word_size)
{
	// Use zlib to read both compressed and uncompressed fasta files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
	
		cerr << "Unable to open file: " << m_filename << endl;
		throw __FILE__ ":parse_fasta: Unable to open fasta file";
	}
	
	const int buffer_len = 2048;
	char buffer[buffer_len];
	
	string defline;
	
	SequenceZ* seq_ptr = NULL;
	size_t seq_len = 0;
	
	Topology local_topo = UNKNOWN_TOPOLOGY;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if(strchr(buffer, '>') != NULL){
			
			if(seq_len > 0){
				
				if(seq_ptr == NULL){
					throw __FILE__ ":parse_fasta: seq_ptr == NULL";
				}
				
				local_topo = get_topology(m_file_info);
								
				// Circularize genomes. The original gottcha_db.pl script
				// added k - 1 bases to the front and k - 1 bases to the back. However,
				// we only need to add k - 1 bases to the back (adding to both back and 
				// front is redundant!).
				if(local_topo == CIRCULAR){
					seq_ptr->circularize(m_word_size);
				}
				
				// The memory associated with this sequence must
				// be deallocated by the calling function
				m_sequence.push_back( make_pair(defline, seq_ptr) );
				
				seq_ptr = NULL;
				seq_len = 0;
			}
			
			// Allocate a new sequence
			seq_ptr = new SequenceZ(m_word_size);

			if(seq_ptr == NULL){
				throw __FILE__ ":parse_fasta: Unable to allocate sequence memory";
			}
			
			// Remove any end of line symbols
			for(char* p = buffer;*p != '\0';++p){
				if( (*p == '\n') || (*p == '\r') ){
					*p = '\0';
				}
			}

			// In cases where the defline is longer than the buffer, indicate truncation
			// with an ellipsis
			if( strlen(buffer) == (buffer_len - 1) ){
				defline = buffer + string("...");
			}
			else{
				defline = buffer;
			}

		}
		else{
			for(char* p = buffer;*p != '\0';++p){
				
				if( !isspace(*p) ){
					
					seq_ptr->push_back(*p);
					++seq_len;
				}
			}
		}
	}
	
	if(seq_len > 0){
		
		if(seq_ptr == NULL){
			throw __FILE__ ":parse_fasta: seq_ptr == NULL";
		}

		local_topo = get_topology(m_file_info);
		
		// Circularize genomes. The original gottcha_db.pl script
		// added k - 1 bases to the front and k - 1 bases to the back. However,
		// we only need to add k - 1 bases to the back (adding to both back and 
		// front is redundant!).
		if(local_topo == CIRCULAR){
			seq_ptr->circularize(m_word_size);
		}

		// The memory associated with this sequence must
		// be deallocated by the calling function
		m_sequence.push_back( make_pair(defline, seq_ptr) );

		seq_ptr = NULL;
		seq_len = 0;
	}
	
	if(seq_ptr != NULL){
		delete seq_ptr;
	}
	
	gzclose(fin);
}

Topology get_topology(const SequenceFileInfo &m_file_info)
{
	if(m_file_info.topo == UNKNOWN_TOPOLOGY){
		throw __FILE__ ":get_topology: Found an unknown sequence topology";
	}
	
	return m_file_info.topo;
}

// Split the defline at the first space after the accession: defline -> {prefix, suffix}
pair<string, string> split_defline(const string &m_defline)
{
	pair<string, string> ret;
	
	const size_t len = m_defline.size();
	
	size_t start = 0;
	
	// Skip any leading space
	while( (start < len) && isspace(m_defline[start]) ){
		++start;
	}
	
	// Make sure that the first non-space character is the fasta header symbol, '>'
	if( (start >= len) || (m_defline[start] != '>') ){
		throw __FILE__ ":split_defline: Malformed fasta header (1)";
	}
	
	size_t stop = start + 1; // Skip the '>'
	
	// Skip any white space between the '>' and the accession
	while( (stop < len) && isspace(m_defline[stop]) ){
		++stop;
	}
	
	if(stop >= len){
		throw __FILE__ ":split_defline: Malformed fasta header (2)";
	}
	
	// The prefix now include all text until the next white space
	while( (stop < len) && !isspace(m_defline[stop]) ){
		++stop;
	}

	// Don't include any trailing '|' in the prefix (since we will be adding a '|'
	// when we modify the defline prefix).
	size_t last = stop;
	
	while(m_defline[last - 1] == '|'){
		--last;
	}
	
	ret.first = m_defline.substr(start, last - start);
	
	// Find the suffix (if any)
	start = stop;
	
	// Skip any white space between the prefix and suffix
	while( (start < len) && isspace(m_defline[start]) ){
		++start;
	}
	
	// The remainder of the defline is the suffix
	ret.second = m_defline.substr(start, len - start);
	
	return ret;
}
