#include <sstream>
#include <iostream>
#include <algorithm>

#include <zlib.h>
#include <string.h>

#include "gottcha_db.h"

using namespace std;

void parse_genbank(const string &m_filename, 
	deque< pair<std::string, SequenceZ*> > &m_sequence, 
	const SequenceFileInfo &m_file_info, const size_t& m_word_size)
{

	// Use zlib to read both compressed and uncompressed genbank files.
	gzFile fin = gzopen(m_filename.c_str(), "r");
	
	if(fin == NULL){
		throw __FILE__ ":parse_genbank: Unable to open genbank file";
	}
		
	const int buffer_len = 2048;
	char buffer[buffer_len];
	char *ptr = NULL;
	const size_t indent = 12;
	
	string definition;
	string accession;
	
	while( gzgets(fin, buffer, buffer_len) ){
		
		if( strlen(buffer) >= (buffer_len -1) ){
			throw __FILE__ ":parse_genbank: Buffer overflow";
		}
		
		// The line that starts with "//" terminates the GBK record
		if(strstr(buffer, "//") == buffer){
			
			// This is an empty record!
			cerr << "Error parsing " << m_filename << endl;
			throw __FILE__ ":parse_genbank: No sequence data found";
		}
		
		if( ( ptr = strstr(buffer, "DEFINITION") ) == buffer){
			
			// Remove any end-of-line symbols
			for(char *i = buffer;*i != '\0';++i){
				
				if( (*i == '\n') || (*i == '\r') ){
					*i = '\0';
				}
			}
			
			ptr += strlen("DEFINITION");
			
			while( (*ptr != '\0') && isspace(*ptr) ){
				++ptr;
			}
			
			if(*ptr == '\0'){
				
				cerr << "Error parsing " << m_filename << endl;
				throw __FILE__ ":parse_genbank: Unable to parse definition";
			}
			
			definition = ptr;
			
			// Handle multi-line definitions
			while( gzgets(fin, buffer, buffer_len) ){
				
				// Remove any end of line symbols
				for(char *i = buffer;*i != '\0';++i){

					if( (*i == '\n') || (*i == '\r') ){
						*i = '\0';
					}
				}

				if(strlen(buffer) <= indent){
					break;
				}
				
				bool valid = true;
				
				for(size_t i = 0;i < indent;++i){
					
					if(buffer[i] != ' '){
						valid = false;
					}
				}
				
				if(!valid){
					break;
				}
				
				definition += (buffer + indent);
			}
			
			
			// Fasta files do not include the period at the end of the definition line
			const size_t def_len = definition.size();
			
			// For compilers that don't support string::back() and string::pop_back()
			if( (def_len > 0) && (definition[def_len - 1] == '.') ){
				definition = definition.substr(0, def_len - 1);
			}
			
			continue;
		}
		
		if( ( ptr = strstr(buffer, "VERSION") ) == buffer){
			
			// Remove any end of line symbols
			for(char *i = buffer;*i != '\0';++i){
				
				if( (*i == '\n') || (*i == '\r') ){
					*i = '\0';
				}
			}
			
			ptr += strlen("VERSION");
			
			// Skip any white space
			while( (*ptr != '\0') && isspace(*ptr) ){
				++ptr;
			}
			
			if(*ptr == '\0'){
			
				cerr << "Error parsing " << m_filename << endl;
				throw __FILE__ ":parse_genbank: Unable to find accession";
			}
			
			while( (*ptr != '\0') && !isspace(*ptr) ){
			
				accession.push_back(*ptr);
				++ptr;
			}

			continue;
		}
		
		if( strstr(buffer, "ORIGIN") == buffer){
			
			SequenceZ* seq_ptr = new SequenceZ(m_word_size);
			
			if(seq_ptr == NULL){
				throw __FILE__ ":parse_genbank: Unable to allocate sequence memory";
			}
			
			size_t seq_len = 0;
	
			// Read the sequence data
			while( gzgets(fin, buffer, buffer_len) ){
				
				if(strstr(buffer, "//") == buffer){
			
					if(seq_len == 0){

						cerr << "Error parsing " << m_filename << endl;
						throw __FILE__ ":parse_genbank: Did not read any sequence data";
					}

					if( accession.empty() ){

						cerr << "Error parsing " << m_filename << endl;
						throw __FILE__ ":parse_genbank: Did not read a valid accession";
					}

					if( definition.empty() ){

						cerr << "Error parsing " << m_filename << endl;
						throw __FILE__ ":parse_genbank: Did not read a valid definition";
					}
					
					const Topology local_topo = get_topology(m_file_info);
										
					// Circularize genomes. The original gottcha_db.pl script
					// added k - 1 bases to the front and k - 1 bases to the back. However,
					// we only need to add k - 1 bases to the back (adding to both back and 
					// front is redundant!).
					if(local_topo == CIRCULAR){
						seq_ptr->circularize(m_word_size);
					}

					stringstream defline;

					// Assume that any accession with an underscore ("_") is a 
					// refseq sequence, gb otherwise
					if(accession.find('_') != string::npos){
						defline << ">ref|" << accession << "| " 
							<< definition;
					}
					else{
						defline << ">gb|" << accession << "| " 
							<< definition;
					}

					m_sequence.push_back( make_pair(defline.str(), seq_ptr) );
					
					seq_ptr = NULL;
					seq_len = 0;
					
					definition.clear();
					accession.clear();

					break;
				}
				
				if(seq_ptr == NULL){
					throw __FILE__ ":parse_genbank: seq_ptr == NULL";
				}
				
				for(ptr = buffer;*ptr != '\0';++ptr){
					
					// Skip digits and white space
					if( isspace(*ptr) || isdigit(*ptr) ){
						continue;
					}
					
					// Make sure that this is a valid character
					if( !isalpha(*ptr) ){
						
						cerr << "Error parsing " << m_filename << endl;
						throw __FILE__ ":parse_genbank: Invalid base";
					}
					
					seq_ptr->push_back(*ptr);
					++seq_len;
				}
			}
			
			// This GenBank record does not have an ORIGIN field (and therefore does not
			// have any sequence).
			if(seq_ptr != NULL){
			
				delete seq_ptr;
				seq_ptr = NULL;
				seq_len = 0;
			}
		}
	}
	
	gzclose(fin);
}

