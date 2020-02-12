#include "gottcha_db.h"
#include "deque_set.h"

#include <iostream>
#include <zlib.h>

using namespace std;

void write_fasta_file(vector<gzFile> &m_fout, const TaxaId &m_taxa_id, 
	const string &m_name, const SequenceZ &m_seq, const unsigned char *m_shared, 
	const unsigned int &m_min_output_frag);
vector<size_t> sort_targets_by_size(const vector<bool> &m_target_set,
        const vector<int> &m_sequence_location,
        const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db);

struct compare_by_group
{
	inline bool operator()(const pair<GroupInfo, deque<Word> > &m_a, const GroupInfo &m_b) const
	{
		return m_a.first < m_b;
	};
};

struct sort_by_file_size
{
	// Sort (file size, file index) in descending order
	inline bool operator()(const pair<size_t, size_t> &m_a, const pair<size_t, size_t> &m_b) const
	{
		if(m_a.first == m_b.first){
			return m_a.second > m_b.second;
		}

		return m_a.first > m_b.first;
	};
};

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

void write_targets(const vector<bool> &m_targets, 
	const deque< pair< GroupInfo, deque<Word> > > &m_target_words,
	const vector<int> &m_sequence_location,
	const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db,
	UpdateInfo &m_progress, const string &m_input_lcp, const Options &m_opt)
{
	// Each rank will be required to write a single file per batch
	// of mpi_numtasks files
	SequenceFile target;
	unsigned char * shared = NULL;
	
	// Ranks will open a separate file for each taxa level
	vector<gzFile> fout(NUM_TAXA_LEVEL, NULL);
	
	// For load balancing, we ideally want all ranks to be writing
	// files that are of equal size. To assist in this goal, we will
	// sort all targets in order of descending file size.
	//
	// Sorting the targets can provide a factor of 2 speed up in run time!
	// For example, when writing 8 taxonomic levels of RefSeq (circa 2018):
	// 	Unsorted target writing: 1.94x10^5 sec
	// 	Sorted target writing: 4.58x10^4 sec
	
	const vector<size_t> target_files = sort_targets_by_size(m_targets,
		m_sequence_location, m_db);

	const size_t num_target = target_files.size();

	for(size_t i = 0;i < num_target;++i){
		
		const int target_rank = i%mpi_numtasks;
		
		SequenceFile local_target;
		
		request_file(target_files[i], m_sequence_location[ target_files[i] ], 
			local_target, m_db);
			
		if(m_opt.verbose >= Options::INFORMATIVE){

			m_progress << "Writing: "
				<< truncate_filename(local_target.name, 75)
				<< " (" << (100.0*i)/num_target << "%)";
			m_progress.flush();
		}

		deque< pair< GroupInfo, deque<Word> > >::const_iterator w_iter =
			lower_bound( m_target_words.begin(), m_target_words.end(), 
				local_target.info.group_id,
				compare_by_group() );

		if( ( w_iter == m_target_words.end() ) || (w_iter->first != local_target.info.group_id) ){
			throw __FILE__ ":write_targets: Unable to lookup target words";
		}
		
		// The total length of all sequence segments in the current target. This is
		// not a cont variable because we will be decrementing it if (and when) the
		// target_len exceeds the 32-bit INT_MAX value that MPI can send in a single
		// call.
		size_t target_len = local_target.size();
	
		unsigned char *local_shared = new unsigned char [target_len];
		
		if(local_shared == NULL){
			throw __FILE__ ":write_targets: Unable to allocate memory for shared sequence buffer";
		}
		
		// By default, each base is unique (i.e. not shared) at each taxonomic
		// level. Since no shared -> NO == 1 for each bit, we need to initialize
		// all bits to 1.
		memset(local_shared, ALL_NO, sizeof(unsigned char)*target_len);
	
		deque<IndexedWord> query;
		
		size_t query_index = 0;
		
		// It is important for computational efficiency to *combine* all of the indexed words
		// from a single target into a single sorted deque. Otherwise, for heavily fragmented
		// genomes, we end up search many small sets of indexed words against the large list of
		// target words (which is very inefficient).
		for(deque< pair<string, SequenceZ*> >::const_iterator j = local_target.seq.begin();
			j != local_target.seq.end();++j){
			
			digest_indexed_words(query, *(j->second), m_opt.word_size, query_index);
			
			query_index += j->second->size();
		}
		
		SORT( query.begin(), query.end() );
		
		deque<Word>::const_iterator iter = w_iter->second.begin();
	
		for(deque<IndexedWord>::const_iterator j = query.begin();j != query.end();++j){

			// Since the words in query are a subset of the words in target word_table, 
			// we can just increment iter until we "catch up" to the word pointed to 
			// by the query iterator. Note that there can be multiple words in query
			// with the same kmer value.
			while(iter->kmer < j->kmer){
				++iter;			
			}

        		// The start of the word is m_word_size bases upstream in the sequence
			unsigned char* shared_end = local_shared + (j->index + 1);

			for(unsigned char* shared_begin = shared_end - m_opt.word_size;
				shared_begin != shared_end;++shared_begin){
				
				*shared_begin &= iter->taxa;
			}
		}
		
		// Reduce using MAX: NOT_TESTED < NO < YES. However, we have examples of genome
		// files whose length exceeds the INT_MAX size that MPI will allow (due to 
		// the 32 bit API limitations of OpenMPI).
		if(mpi_rank == target_rank){

			unsigned char *ptr = local_shared;
			
			// Reduce in multiple steps to avoid the MPI 32-bit size limitation
			while(target_len > 0){
			
				const size_t chunk_size = min( size_t(INT_MAX >> 1), target_len );
				
				// Reduce with the binary AND operation
				MPI_Reduce(MPI_IN_PLACE, ptr, chunk_size, 
					MPI_BYTE, MPI_BAND, target_rank, MPI_COMM_WORLD);
				
				target_len -= chunk_size;
				ptr += chunk_size;
			}
			
			// The target rank saves a pointer to this buffer (which will be deleted
			// later).
			shared = local_shared;
		}
		else{ // mpi_rank != target_rank

			unsigned char *ptr = local_shared;
			
			// Reduce in multiple steps to avoid the MPI 32-bit size limitation
			while(target_len > 0){
				
				const size_t chunk_size = min( size_t(INT_MAX >> 1), target_len );
				
				// Reduce with the binary AND operation
				MPI_Reduce(ptr, NULL, chunk_size,
					 MPI_BYTE, MPI_BAND, target_rank, MPI_COMM_WORLD);
					
				target_len -= chunk_size;
				ptr += chunk_size;
			}
			
			// We must manually deallocate the shared mask *only* for the
			// *non* target rank
			delete [] local_shared;
			local_shared = NULL;
		}
				
		if(target_rank == mpi_rank){
			
			// The target rank saves a copy of this file
			target = local_target;
		}
		else{
			// Only the non-target rank frees the memory associated with the local target
			// due to shallow copy semantics of SequenceFile;
			local_target.clear();
		}
		
		// Have we processed mpi_numtasks number of sequence files or are we at the last
		// sequence file? If so, it is time to flush to disk.
		if( ( ( target_rank == (mpi_numtasks - 1) ) || (i == (num_target - 1) ) ) ){
		
			if(m_opt.verbose >= Options::INFORMATIVE){

				m_progress << "Flushing files to disk ("
					<< (100.0*i)/num_target << "%)";
				m_progress.flush();
			}

			if( shared != NULL){
			
				#pragma omp parallel for
				for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

					if(target.info.group_id[level] == UNDEFINED_GROUP){
						continue;
					}

					unordered_map<TaxaLevel, string>::const_iterator prefix_iter = 
						m_opt.taxa_level_output_prefix.find( TaxaLevel(level) );

					if( prefix_iter == m_opt.taxa_level_output_prefix.end() ){
						throw __FILE__ ":write_targets: Unable to lookup output file prefix";
					}

					const string filename = format_output_filename(target.name, 
							prefix_iter->second, m_opt.output_root_dir, m_input_lcp,
							m_opt.compress_output);

					// If we are re-rooting the output files, we will likely need to create new directories
					// to ensure that we are writing to a valid path.
					// ** Update ** this directory test is created by the rank 0 process at the
					// start of the calculation.

					fout[level] = gzopen( filename.c_str(), m_opt.compress_output ? "w" : "wT");

					if(fout[level] == NULL){

						cerr << "Unable to open " << filename << " for writing fasta output at level: " 
							<< level << endl;
						throw __FILE__ ":write_targets: Unable to open file";
					}
				}

				size_t sequence_offset = 0;
				
				for(deque< pair<string, SequenceZ*> >::const_iterator j = target.seq.begin();
					j != target.seq.end();++j){

					// Each rank writes one or more fragmented fasta files. Note that we need to pass
					// the taxa id to be included in the fasta defline (so that the GOTTCHA profile
					// tool has access to the taxa id for each sequence).
					write_fasta_file(fout, target.info.taxa_id, 
						j->first, *(j->second), 
						shared + sequence_offset,
						m_opt.min_output_frag);
						
					sequence_offset += j->second->size();
				}

				delete [] shared;
				shared = NULL;
				
				target.clear();
				
				#pragma omp parallel for
				for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

					if(fout[level] != NULL){
					
						gzclose(fout[level]);
						fout[level] = NULL;
					}
				}
			}
			
			// Use an explicit barrier here to synchronize all ranks. Since some ranks may have
                        // much more work to do (i.e. writing a large fasta file to disk), the other ranks can
                        // get ahead if the MPI library starts caching the MPI_Reduce
                        // results. This could result is extra memory consumption that exhausts the
                        // memory on the ranks that are writing large files.
                        MPI_Barrier(MPI_COMM_WORLD);
		}
	}
}

void write_fasta_file(vector<gzFile> &m_fout, const TaxaId &m_taxa_id, 
	const string &m_name, const SequenceZ &m_seq, 
	const unsigned char *m_shared, const unsigned int &m_min_output_frag)
{
	const size_t fasta_chunk = 70; // Write fasta sequences with 70 columns per line
	
	// Split the defline at the first space after the accession
	const pair<string, string> defline = split_defline(m_name);
		
	#pragma omp parallel for
	for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

		gzFile &fout = m_fout[level];

		if(fout == NULL){
			continue;
		}
				
		deque<char> fragment;

		bool valid_fragment = false;
		unsigned int start = 0;
		unsigned int stop = 0;

		unsigned int index = 0;
		SequenceZ::const_iterator i = m_seq.begin();

		while(true){

			bool end_of_seq = ( i == m_seq.end() );

			const char base = (end_of_seq ? '?' : *i);
			const unsigned char shared = (end_of_seq ? ALL_NO : m_shared[index]);

			if( !end_of_seq && ( ( (shared >> level) & 1 ) == NO ) ){ // Not shared

				// This is a signature base! Append to the current fragment, or start
				// a new fragment
				fragment.push_back(base);

				if(valid_fragment){
					++stop;
				}
				else{
					valid_fragment = true;

					// For the acutal start and stop coordinates, we need to
					// query the SequenceZ iterator (since non-ATGC bases are
					// skipped in the stored sequence, which causes offsets in
					// the base counting).
					start = stop = i.index();
				}
			}
			else{ // A shared (i.e. non-signature) base, or the end of the sequence

				// Since we have found an invalid base, terminate any existing fragment
				if(valid_fragment){

					valid_fragment = false;

					const unsigned int frag_len = stop - start + 1;

					if(frag_len >= m_min_output_frag){

						// The fragment coordinates are 1-based
						// Please note that gzprintf can only write a maximum of 8191 bytes
						// of uncompressed data!

						// Write the fasta defline
						gzputs( fout, defline.first.c_str() );
						gzprintf(fout, "|%d|%d|%s| ",
							start + 1, stop + 1, 
							m_taxa_id.c_str() );
						gzputs( fout, defline.second.c_str() );
						gzputc(fout, '\n');

						// Write the fasta sequence
						deque<char>::const_iterator frag_iter = fragment.begin();
						
						for(size_t j = 0;j < frag_len;j += fasta_chunk){

							const size_t local_len = min(frag_len - j, fasta_chunk);

							for(size_t k = 0;k < local_len;++k, ++frag_iter){
								gzputc(fout, *frag_iter);
							}

							gzputc(fout, '\n');
						}
					}

					fragment.clear();
				}
			}

			if(end_of_seq){
				break;
			}

			// Increment *after* we test for the end of sequence (since incrementing
			// past end() will generate an error).
			++i;
			++index;
		}
	}
}

// To assist in load balancing, we will be writing files in order of 
// descending file size.
vector<size_t> sort_targets_by_size(const vector<bool> &m_target_set,
	const vector<int> &m_sequence_location, 
	const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db)
{
	const size_t num_file = m_target_set.size();
	
	size_t num_valid_file = 0;

	// Count the number of valid files in the target set. We need this
	// count to allocate a buffer for MPI AllReduce (we could have used
	// MPI Bcast as a quick and dirty solution, but this would 
	// potentially required many calls to Bcast).
	for(size_t i = 0;i < num_file;++i){
		num_valid_file += m_target_set[i] ? 1 : 0;
	}

	unsigned long long int *buffer = new unsigned long long int [num_valid_file];

	if(buffer == NULL){
		throw __FILE__ ":sort_targets_by_size: Unable to allocate buffer";
	}

	size_t index = 0;

	for(size_t i = 0;i < num_file;++i){

		if(m_target_set[i]){

			buffer[index] = 0;
	
			if(m_sequence_location[i] == mpi_rank){

				unordered_map<size_t, pair<unsigned char*, size_t> >::const_iterator iter =
					m_db.find(i);

				if( iter == m_db.end() ){
					throw __FILE__ ":sort_targets_by_size: Unable to find sequence";
				}

				buffer[index] = iter->second.second;
			}

			++index;
		}
	}

	MPI_Allreduce(MPI_IN_PLACE, buffer, num_valid_file, 
		MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD);

	vector< pair<size_t, size_t> > file_size(num_valid_file);

	index = 0;

	for(size_t i = 0;i < num_file;++i){
                
                if(m_target_set[i]){
		
			// Recall that: 
			//	"index" is relative to the buffer
			//	"i" is the actual file id
			file_size[index] = make_pair(buffer[index], i);
			++index;
		}
	}

	delete [] buffer;

	// All ranks sort the file size data on their own. Note that the
	// comparison function sorts on file size as the primary key and
	// file index as the secondary key, so the sorted order is
	// uniquely defined so all ranks will agree on the order of
	// files.
	SORT( file_size.begin(), file_size.end(), sort_by_file_size()  );

	// After sorting by file size in descending order. we no longer need
	// to store the file size (only return the file index).	
	vector<size_t> ret(num_valid_file);

	for(size_t i = 0;i < num_valid_file;++i){
		ret[i] = file_size[i].second;
	}

	return ret;
}
