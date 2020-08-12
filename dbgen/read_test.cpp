// Test bed for I/O benchmarking
// J. D. Gans
// Bioscience Division, B-11
// Los Alamos National Laboratory
// Wed Mar 21 10:44:05 2018

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include <stdlib.h>
#include <stdlib.h>
#include <mpi.h>
#include <signal.h>
#include <math.h>
#include <string.h> // memcpy

using namespace std;

// Global variables
int mpi_rank;
int mpi_numtasks;
double start_time;

void terminate_program(int m_sig);
string report_run_time();
deque<string> tab_split(const string &m_path);
deque<string> parse_fasta(const string &m_filename);

int main(int argc, char *argv[])
{
	try{
		
		MPI_Init(&argc, &argv);
		MPI_Comm_size(MPI_COMM_WORLD, &mpi_numtasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
		
		signal( SIGINT, terminate_program );
		signal( SIGTERM, terminate_program );
		signal( SIGSEGV, terminate_program );
		
		start_time = MPI_Wtime();
		
		const string root = "/users/218817/scratch/projects/20171102_gottcha2_db/genomes";
		
		cerr << "I/O testbed [" << mpi_rank << "]" << endl;
		
		// The filename of target files
		const string genome_mapping_file = 
			"/users/218817/scratch/projects/20171102_gottcha2_db/genomes/RefSeq_genomes.strain.list";
		
		deque<string> target_files;
		
		ifstream fin( genome_mapping_file.c_str() );
		
		if(!fin){
			throw __FILE__ ":main: Unable to open genome mapping file for reading";
		}
		
		string line;
		
		double profile = MPI_Wtime();
		
		while( getline(fin, line) ){
			
			const deque<string> data = tab_split(line);
			
			if(data.size() > 2){
				target_files.push_back( data[1] );
			}
			
			if(target_files.size() > 10000){
				break;
			}
			
		}
		
		profile = MPI_Wtime() - profile;
		
		cerr << "[" << mpi_rank << "] Found " << target_files.size() << " target files in " 
			<< profile << " sec" << endl;
		
		size_t num_seq = 0;
		size_t num_base = 0;
		
		profile = MPI_Wtime();
		
		//#define	PARALLEL_READ
		
		#ifdef PARALLEL_READ
		for(deque<string>::const_iterator i = target_files.begin();i != target_files.end();++i){
			
			const deque<string> seq = parse_fasta(root + "/" + *i);
			
			num_seq += seq.size();
			
			for(deque<string>::const_iterator j = seq.begin();j != seq.end();++j){
				num_base += j->size();
			}
		}
		#else // #ifdef PARALLEL_READ --> SERIAL_READ
		for(deque<string>::const_iterator i = target_files.begin();i != target_files.end();++i){
			
			deque<string> seq;
			
			unsigned int buffer_size = 0;
			unsigned char* buffer = NULL;
			
			if(mpi_rank == 0){
				
				seq = parse_fasta(root + "/" + *i);
				
				// Send this sequence data to the workers
				buffer_size = sizeof(unsigned int);
				
				for(deque<string>::const_iterator j = seq.begin();j != seq.end();++j){
					buffer_size += j->size() + 1; // Include the '\0'
				}
				
				MPI_Bcast( (void*)&buffer_size, sizeof(unsigned int), MPI_BYTE, 0, MPI_COMM_WORLD );
				
				buffer = new unsigned char [buffer_size];
				
				if(buffer == NULL){
					throw __FILE__ ":main: Unable to allocate sequence buffer [0]";
				}
				
				unsigned char* ptr = buffer;
				
				unsigned int N = seq.size();
				
				memcpy( ptr, &N, sizeof(unsigned int) );
				ptr += sizeof(unsigned int);
				
				for(deque<string>::const_iterator j = seq.begin();j != seq.end();++j){
					
					memcpy( ptr, j->c_str(), j->size() + 1 );
					ptr += j->size() + 1;
				}
				
				MPI_Bcast( (void*)buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD );
				
				delete [] buffer;
			}
			else{
				
				MPI_Bcast( (void*)&buffer_size, sizeof(unsigned int), MPI_BYTE, 0, MPI_COMM_WORLD );
				
				buffer = new unsigned char [buffer_size];
				
				if(buffer == NULL){
					throw __FILE__ ":main: Unable to allocate sequence buffer [!0]";
				}
				
				MPI_Bcast( (void*)buffer, buffer_size, MPI_BYTE, 0, MPI_COMM_WORLD );
				
				unsigned char* ptr = buffer;
				
				unsigned int N = 0;
				
				memcpy( &N, ptr, sizeof(unsigned int) );
				ptr += sizeof(unsigned int);
				
				for(unsigned int j = 0;j < N;++j){
					
					seq.push_back( (char*)ptr );
					ptr += seq.back().size() + 1;
					
				}
				
				delete [] buffer;
			}
			
			num_seq += seq.size();
			
			for(deque<string>::const_iterator j = seq.begin();j != seq.end();++j){
				num_base += j->size();
			}
		}
		
		#endif // PARALLEL_READ
		
		profile = MPI_Wtime() - profile;
		
		cerr << "[" << mpi_rank << "] read " << num_base 
			<< " bases in " << num_seq << " seqs in " << profile << " sec" << endl;
		
		MPI_Finalize();
		
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

deque<string> parse_fasta(const string &m_filename)
{
	ifstream fin( m_filename.c_str() );
	
	if(!fin){
		cerr << "[" << mpi_rank << "] Unable to open: " << m_filename << endl;
		throw __FILE__ ":parse_fasta: Unable to open fasta file";
	}
	
	deque<string> ret;
	
	string line;
	string seq;
	
	while( getline(fin, line) ){
		
		// Skip the fasta headers
		if( line.find('>') != string::npos ){
		
			if( !seq.empty() ){
			
				ret.push_back(seq);
				seq.clear();
			}
		}
		else{
			
			seq += line;
		}
	}
	
	if( !seq.empty() ){
			
		ret.push_back(seq);
		seq.clear();
	}
			
	return ret;
}

deque<string> tab_split(const string &m_path)
{
	const char m_delim = '\t';
	
	deque<string> ret;
	
	string curr;
	
	for(string::const_iterator i = m_path.begin();i != m_path.end();++i){
		
		if(*i == m_delim){
			
			if( !curr.empty() ){
				ret.push_back(curr);
			}
			
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

// Report our rank, the signal we caught and the time spent running
void terminate_program(int m_sig)
{
	cerr << "[" << mpi_rank << "] caught signal " << m_sig << endl;
	cerr << report_run_time() << endl;
	
	MPI_Abort(MPI_COMM_WORLD, 0);
}

// Run time computes the total run time. The results are formatted as a string.
string report_run_time()
{
	double elapsed_time = MPI_Wtime() - start_time; // In sec
	
	const double elapsed_sec = fmod(elapsed_time, 60.0);
	
	elapsed_time = (elapsed_time - elapsed_sec)/60.0; // In min
	
	const double elapsed_min = fmod(elapsed_time, 60.0);
	elapsed_time = (elapsed_time - elapsed_min)/60.0; // In hour
	
	const double elapsed_hour = fmod(elapsed_time, 24.0);
	elapsed_time = (elapsed_time - elapsed_hour)/24.0; // In day
	
	stringstream sout;
	
	sout << "Run time is " 
		<< elapsed_time 
		<< " days, "
		<< elapsed_hour
		<< " hours, "
		<< elapsed_min
		<< " min and "
		<< elapsed_sec
		<< " sec";
	
	return sout.str();
}
