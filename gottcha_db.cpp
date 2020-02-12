// gottcha_db
// A parallel version of the original (serial, perl) gottcha_db program
// J. D. Gans
// Bioscience Division, B-10
// Los Alamos National Laboratory
// June, 2018

#include <iostream>
#include <sstream>
#include <iomanip>

#include <stdlib.h>
#include <mpi.h>
#include <signal.h>
#include <math.h>

#include "gottcha_db.h"
#include "deque_set.h"

using namespace std;

// Global variables for MPI
int mpi_numtasks;
int mpi_rank;

double start_time;

void terminate_program(int m_sig);
string report_run_time();

size_t per_rank_database_size(const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db);
	
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
				
		// Only rank 0 loads the command line options
		Options opt(argc, argv, mpi_rank == 0);
		
		// Rank 0 shares the commnad line options with the other ranks
		broadcast(opt, mpi_rank /* my rank */, 0 /* src rank*/);
		
		if(opt.quit){
		
			MPI_Finalize();

			return EXIT_SUCCESS;
		}
		
		UpdateInfo progress;
		
		// Only rank 0 keeps the user updated
		if(mpi_rank == 0){
			
			progress.Log(opt.log_file);
		
			// Reduce the log file clutter by only printing 3 significant
			// figures
			progress << setprecision(3);
			
			if(opt.verbose >= Options::INFORMATIVE){

				cerr << "GOTTCHA database creation tool version " << VERSION << endl;
				cerr << "\tRunning on " << mpi_numtasks << " MPI ranks" << endl;
				cerr << "\tWord size = " << opt.word_size << endl;
				cerr << "\tMinimum output fragment size = " << opt.min_output_frag << endl;
				
				if(opt.max_worker_RAM != 0){
					cerr << "\tMemory is limited to " << float(opt.max_worker_RAM)/GB
						<< " GB of RAM" << endl;
				}
				else{
					cerr << "\tMemory is limited to " << RAM_FRACTION*100.0
						<< "% of available RAM" << endl;
				}
				
				if(opt.compress_output){
					cerr << "\tWriting compressed fasta files" << endl;
				}
				else{
					cerr << "\tWriting uncompressed fasta files" << endl;
				}
		
				cerr << endl;
			}
		}
		else{ // mpi_rank > 0
			opt.verbose = Options::SILENT;
		}
		
		// "opt.max_worker_RAM == 0" is a special condition that forcs the compute
		// nodes to measure the amount of available RAM and restrict their memory
		// usage to RAM_FRACTION*available_RAM().
		if(opt.max_worker_RAM == 0){
		
			opt.max_worker_RAM = available_RAM()*RAM_FRACTION;
			
			if(opt.max_worker_RAM == 0){
				throw __FILE__ ":compute_main: No available RAM!";
			}
		}
		
		deque<unsigned int> groups_to_squash;
		
		// On which rank is each sequence file stored?
		vector<int> sequence_location;
		
		/////////////////////////////////////////////////////////////////////////////
		// Load and distribute all sequences into the memory of all ranks
		/////////////////////////////////////////////////////////////////////////////
		
		// The database of SequenceFile's that have been packed into binary buffers for
		// easy MPI transport.
		unordered_map<size_t, pair<unsigned char*, size_t> > db;
		
		// The input_lcp is the longest common prefix of the input files. Only needed by
		// the root when the user has specified a new output root directory
		const string input_lcp = 
			load_sequence_database(db, sequence_location, groups_to_squash, progress, opt);
		
		const size_t num_seq_file = sequence_location.size();
		
		// The amount of RAM consumed by the local database of sequences
		// (which we need to known to make sure an individual rank does
		// not consume too much memory).
		const size_t local_database_size = per_rank_database_size(db);
		
		// Profile the various execution times
		double profile_target_load = 0.0;
		double profile_background_request = 0.0;
		double profile_background_subtract = 0.0;
		double profile_target_write = 0.0;
		
		// Process all of the sequence files
		size_t current_target = 0;

		while(current_target < num_seq_file){
			
			double profile_start = MPI_Wtime();
			
			// We will load target sequences until the first rank runs out of available
			// working RAM
			
			// The active set of target words that have been retrieved from 
			// the distributed database of sequences.
			deque< pair< GroupInfo, deque<Word> > > target_words;
			
			// Track the files that have been loaded into the current target set.
			// These files will *not* be loaded again as background files
			vector<bool> valid_target_set(num_seq_file, false);
			
			// The amount of RAM used to store the active target sequences and associated words
			// on this rank
			size_t target_ram = local_database_size;

			// Use a do-while loop to ensure that at least one target is loaded
			do{
				// Every rank loads the same sequence from the rank that is storing it
				size_t local_size = 0;
				
				// Use a while loop to skip squashed sequences
				while( (local_size == 0) && (current_target < num_seq_file) ){
				
					string target_name;
					
					// Each rank stores a non-overlapping subset of words
					local_size = request_file_words(current_target, 
						sequence_location[current_target],
						target_name, target_words,
						groups_to_squash, db, opt.word_size);
					
					if(local_size > 0){
						
						// If local_size != 0, we have a valid target
						valid_target_set[current_target] = true;
						
						if(opt.verbose >= Options::INFORMATIVE){

							progress << "Target: "
								<< truncate_filename(target_name, 75) 
								<< ": " << (current_target*100.0)/num_seq_file
								<< '%';
							progress.flush();
						}
					}
					
					++current_target;
				}
				
				target_ram += local_size;
			}
			while( (current_target < num_seq_file) && 
				room_for_more(opt.max_worker_RAM > 2*target_ram) ); // Leave room for the background
										    // by doubling the target_ram
										    // when testing for memory
										    // consumption.
			
			size_t num_target_words = 0;
			
			for(deque< pair< GroupInfo, deque<Word> > >::const_iterator i = target_words.begin();
				i != target_words.end();++i){
				
				num_target_words += i->second.size();
			}

			size_t min_words = 0;
			size_t max_words = 0;
			
			MPI_Reduce(&num_target_words, &min_words, 1, 
					MPI_UNSIGNED_LONG_LONG, MPI_MIN, 0, MPI_COMM_WORLD);
			
			MPI_Reduce(&num_target_words, &max_words, 1, 
					MPI_UNSIGNED_LONG_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
			
			const double target_load = MPI_Wtime() - profile_start;
			
			profile_target_load += target_load;
			
			if(opt.verbose >= Options::INFORMATIVE){
				
				const double word_ratio = double(min_words)/max_words;
				
				//progress << "Loaded " << (current_target*100.0)/num_seq_file 
				//	<< "% of targets in " << target_load << " sec (balance = "
				//	<< word_ratio << "): " << num_target_words << " words";
				progress << "Loaded " << (current_target*100.0)/num_seq_file 
					<< "% of targets in " << target_load << " sec (balance = "
					<< word_ratio << ")";
				progress.flush();
			}

			const double target_progress = (current_target*100.0)/num_seq_file;
			
			/////////////////////////////////////////////////////////////
			// Subtract all of the background sequences from the targets
			/////////////////////////////////////////////////////////////
			
			size_t i = 0;
			
			// The first time subtracting sequences, we will need to subtract
			// targets from targets
			bool subtract_target = true;
			
			while(i < num_seq_file){
			
				size_t background_ram = target_ram;
				
				deque< pair< GroupInfo, deque<Word> > > background_words;

				profile_start = MPI_Wtime();

				// The following while loop accumulates background sequences with identical
				// taxonomic groupings (i.e. the group_id field in SequenceFileInfo is the
				// same).
				size_t num_background = 0;

				// Use a do-while loop to ensure that at least one background is loaded
				do{
					// Don't reload files that are in the current valid target
					// set, as these have already been subtracted from
					// each other
					if(valid_target_set[i]){

						// Skip this sequence file, it is in the
						// current target set
						++i;
						continue;
					}

					++num_background;

					string background_name;

					// Each rank stores a non-overlapping subset of words
					background_ram += request_file_words(i, 
						sequence_location[i],
						background_name, background_words,
						deque<unsigned int>(), db, opt.word_size);

					++i;
					
					if(opt.verbose >= Options::VERBOSE){

						progress << "Background: "
							<< truncate_filename(background_name, 75) 
							<< ": " << (i*100.0)/num_seq_file
							<< '%';
						progress.flush();
					}
				}
				while( (i < num_seq_file) && 
					room_for_more(opt.max_worker_RAM > background_ram) );

				size_t num_background_words = 0;

				for(deque< pair< GroupInfo, deque<Word> > >::const_iterator j = 
					background_words.begin();j != background_words.end();++j){

					num_background_words += j->second.size();
				}

				const double background_request = MPI_Wtime() - profile_start;

				profile_background_request += background_request;

				if( (num_background > 0) && (opt.verbose >= Options::INFORMATIVE) ){

					//progress << target_progress << "%: Loaded " 
					//	<< num_background << " background in " 
					//	<< background_request << " sec: "
					//	<< num_background_words << " words";
					progress << target_progress << "%: Loaded " 
						<< num_background << " background in " 
						<< background_request << " sec";

					progress.flush();
				}

				profile_start = MPI_Wtime();

				// Each rank subtracts a non-overlapping subset of words
				// Subtract the background words from all of the target sequences
				// and, if subtract_target == true, subtract the target sequences
				// from themselves
				subtract_sequence(target_words, background_words, 
					progress, (opt.verbose >= Options::INFORMATIVE),
					subtract_target);
				
				// Target vs target subtraction only needs to happen once for each 
				// set of targets
				subtract_target = false;
				
				const double background_subtract = MPI_Wtime() - profile_start;

				profile_background_subtract += background_subtract;

				if( (num_background > 0) && (opt.verbose >= Options::INFORMATIVE) ){

					//progress << target_progress << "%: Background X " 
					//	<< num_background << " in " 
					//	<< background_subtract << " sec: "
					//	<< num_background_words << " words (" 
					//	<< (num_target_words + num_background_words)/background_subtract 
					//	<< " w/s): " 
					//	<< (i*100.0)/num_seq_file
					//	<< '%';
					progress << target_progress << "%: Subtract background X " 
						<< num_background << " in " 
						<< background_subtract << " sec: "
						<< (i*100.0)/num_seq_file
						<< '%';
						
					progress.flush();
				}
			}

			profile_start = MPI_Wtime();
			
			// Write the subtracted target sequences (which will require collecting
			// the unique signature information from all ranks, since they all store
			// non-overlapping subsets of words)
			
			write_targets(valid_target_set, target_words, 
				sequence_location, db,
				progress, input_lcp, opt);
			
			profile_target_write += MPI_Wtime() - profile_start;
		}
		
		if(opt.verbose >= Options::INFORMATIVE){

			progress << "Loaded targets in " << profile_target_load << " sec" << endl;
			progress << "Requested background in " << profile_background_request << " sec" << endl;
			progress << "Subtracted background & targets in " << profile_background_subtract << " sec" << endl;
			progress << "Wrote targets in " << profile_target_write << " sec" << endl;
			progress.flush();
		}

		// We need to manually clean up all of the allocated sequence data buffers
		for(unordered_map< size_t, pair<unsigned char*, size_t> >::iterator i = db.begin();i != db.end();++i){
			delete [] i->second.first;
		}
		
		MPI_Finalize();
		
		if(mpi_rank == 0){
			cerr << report_run_time() << endl;
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

