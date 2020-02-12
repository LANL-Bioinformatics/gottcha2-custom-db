#include "gottcha_db.h"
#include "deque_set.h"

#include <omp.h>

using namespace std;

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

void digest_words(deque<Word> &m_words, const SequenceZ &m_seq,
	const unsigned int &m_word_size, 
	const TaxaSet &m_init_taxa /* = NO_TAXA*/)
{
	const BaseWord comp_shift = 2*(m_word_size - 1);
	const BaseWord mask = ( 1UL << (2*m_word_size) ) - 1;

	const size_t len = m_seq.size();
	
	#pragma omp parallel
	{
		const size_t num_thread = omp_get_num_threads();
		const size_t tid = omp_get_thread_num();

		const size_t start = (tid*len)/num_thread;
		const size_t stop = min( ( (tid + 1)*len )/num_thread + m_word_size, len);
		
		const SequenceZ::const_iterator seq_begin = m_seq.begin() + start;
		const SequenceZ::const_iterator seq_end = m_seq.begin() + stop;
		
		BaseWord w = 0;
		BaseWord comp = 0;
		unsigned int word_len = 0;

		Word address(0, m_init_taxa);
		deque<Word> local;
		
		for(SequenceZ::const_iterator i = seq_begin;i != seq_end;++i){

			// Restart the word every time we start a new fragment. Otherwise,
			// simply increment the word_len.
			word_len = i.new_fragment() ? 1 : word_len + 1;

			const unsigned char base = i.binary_base();

			w = (w << 2) | base;

			// We have defined the binary values of A, T, G and C such that
			// T - base = complement(base)

			comp = (comp >> 2) | (BaseWord(SequenceZ::T - base) << comp_shift);

			address.kmer = min(w & mask, comp & mask);

			// Don't use "min(w & mask, comp & mask)" to compute the destination
			// when distributing words to ranks or threads
			// (it leads to an unequal distribution -- some partitions getting
			// a lot of words to process, some getting just a few). Using the
			// exclusive or (^) appears to give a much more even distribution of
			// words (which will lead to better load balancing and more efficient
			// use of memory).
			//Word address = ( (w ^ comp) & mask ) % mpi_numtasks;
			//if( ( (w ^ comp) & mask ) % mpi_numtasks == Word(mpi_rank) ){
			//
			// *update*: It turns out that designing a really good hash function is an
			// art form. See the kmer_hash() for more information!
			if( (word_len >= m_word_size) && 
				(kmer_hash(address.kmer & RANK_HASH_BITS, mpi_numtasks) == mpi_rank) ){

				local.push_back(address);
			}
		}
		
		#pragma omp critical
		for(deque<Word>::const_iterator i = local.begin();i != local.end();++i){
			m_words.push_back(*i);
		}
	}	
}

void digest_indexed_words(deque<IndexedWord> &m_words, const SequenceZ &m_seq, 
	const unsigned int &m_word_size, const size_t &m_index)
{
	const BaseWord comp_shift = 2*(m_word_size - 1);
	const BaseWord mask = ( 1UL << (2*m_word_size) ) - 1;

	const size_t len = m_seq.size();
	
	#pragma omp parallel
	{
		const size_t num_thread = omp_get_num_threads();
		const size_t tid = omp_get_thread_num();
		
		const size_t start = (tid*len)/num_thread;
		const size_t stop = min( ( (tid + 1)*len )/num_thread + m_word_size, len);
		
		const SequenceZ::const_iterator seq_begin = m_seq.begin() + start;
		const SequenceZ::const_iterator seq_end = m_seq.begin() + stop;
		
		BaseWord w = 0;
		BaseWord comp = 0;
		unsigned int word_len = 0;

		IndexedWord address(0, m_index + start);
		deque<IndexedWord> local;
		
		for(SequenceZ::const_iterator i = seq_begin;i != seq_end;++i, ++address.index){

			// Restart the word every time we start a new fragment. Otherwise,
			// simply increment the word_len.
			word_len = i.new_fragment() ? 1 : word_len + 1;

			const unsigned char base = i.binary_base();

			w = (w << 2) | base;

			// We have defined the binary values of A, T, G and C such that
			// T - base = complement(base)

			comp = (comp >> 2) | (BaseWord(SequenceZ::T - base) << comp_shift);

			address.kmer = min(w & mask, comp & mask);

			// Don't use "min(w & mask, comp & mask)" to compute the destination
			// when distributing words to ranks or threads
			// (it leads to an unequal distribution -- some partitions getting
			// a lot of words to process, some getting just a few). Using the
			// exclusive or (^) appears to give a much more even distribution of
			// words (which will lead to better load balancing and more efficient
			// use of memory).
			//Word address = ( (w ^ comp) & mask ) % mpi_numtasks;
			//if( ( (w ^ comp) & mask ) % mpi_numtasks == Word(mpi_rank) ){
			//
			// *update*: It turns out that designing a really good hash function is an
			// art form. See the kmer_hash() for more information!
			if( (word_len >= m_word_size) && 
				(kmer_hash(address.kmer & RANK_HASH_BITS, mpi_numtasks) == mpi_rank) ){

				local.push_back(address);
			}
		}
		
		#pragma omp critical
		for(deque<IndexedWord>::const_iterator i = local.begin();i != local.end();++i){
			m_words.push_back(*i);
		}
	}
}

size_t request_file_words(const size_t &m_target_index, const int &m_target_rank,
	string &m_target_name, 
	deque< pair< GroupInfo, deque<Word> > > &m_target_words,
	const deque<unsigned int> &m_groups_to_squash, 
	const unordered_map<size_t, pair<unsigned char*, size_t> > &m_db,
	const size_t &m_word_size)
{
	unsigned char *buffer = NULL;
	size_t buffer_size = 0;
	TaxaSet init_taxa;
	
	// Every rank needs to load the current target which is stored by m_target_rank
	if(mpi_rank == m_target_rank){
	
		unordered_map<size_t, pair<unsigned char*, size_t> >::const_iterator iter = 
			m_db.find(m_target_index);
		
		if( iter == m_db.end() ){
			throw __FILE__ ":request_file_words: Unable to find requested sequence file";
		}
		
		buffer = iter->second.first;
		buffer_size = iter->second.second;
		
		// Is this sequence entirely squashed?
		SequenceFileInfo info;
			
		mpi_unpack(buffer, info);
			
		init_taxa = get_taxa_mask(m_groups_to_squash, info);
		
		// If the entire file has been squashed, set the buffer size to 0 so that the
		// other ranks will not expect to load sequence data
		if(init_taxa == NO_TAXA){
			buffer_size = 0;
		}
		
		MPI_Bcast((void*)&buffer_size, sizeof(size_t), MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
		
		if(buffer_size == 0){
			return 0;
		}
		
		MPI_Bcast(buffer, iter->second.second, MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
		MPI_Bcast(&init_taxa, sizeof(TaxaSet), MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
	}
	else{
		
		MPI_Bcast(&buffer_size, sizeof(size_t), MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
		
		if(buffer_size == 0){
			return 0;
		}
		
		buffer = new unsigned char [buffer_size];
		
		if(buffer == NULL){
			throw __FILE__ ":request_file_words: Unable to allocate buffer";
		}
		
		MPI_Bcast(buffer, buffer_size, MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
		MPI_Bcast(&init_taxa, sizeof(TaxaSet), MPI_BYTE, m_target_rank, MPI_COMM_WORLD);
		
	}
	
	SequenceFile local;
	
	mpi_unpack(buffer, local);
	
	if(mpi_rank != m_target_rank){
		delete [] buffer;
	}
	
	m_target_name = local.name;
	
	const size_t num_seq = local.seq.size();

	// Assume that the target files are loaded in order of ascending group_id
	if( m_target_words.empty() || (m_target_words.back().first != local.info.group_id) ){
		m_target_words.push_back( make_pair( local.info.group_id, deque<Word>() ) );
	}
	
	deque<Word> &w = m_target_words.back().second;
	
	deque<Word> local_w;
	
	// Accumulate the words for this target into a local data structure to
	// avoid the cost of resorting words that were sorted on previous iterations
	for(size_t i = 0;i < num_seq;++i){
		
		digest_words(local_w, *(local.seq[i].second), m_word_size, init_taxa);
		
		// Free the sequence data as we go
		delete local.seq[i].second;
	}
	
	// The default parallel sort algorithm used in __gnu_parallel::sort (potentially 
	// called via the SORT macro) is the "multiway merge sort" (MWMS) and 
	// requires *twice* the memory of input array to sort in parallel! The other 
	// sort algorithm implemented by __gnu_parallel::sort, i.e. quick sort 
	// (and balanced quick sort), does not require any extra RAM, but offers 
	// much worse parallel speed up (two fold at best, unless nested parallelism
	// is enabled for OpenMP).
	// * Force MWMS using: __gnu_parallel::multiway_mergesort_tag()
	// * Force balanced QS using: __gnu_parallel::balanced_quicksort_tag()

	SORT( local_w.begin(), local_w.end(), __gnu_parallel::balanced_quicksort_tag() );
	
	local_w.erase( unique( local_w.begin(), local_w.end() ), local_w.end() );
	
	const size_t init_word_count = w.size();
	
	if( w.empty() ){
		
		// The first target that matches this taxonomic grouping we can
		// simply copy the words
		w.swap(local_w);
	}
	else{
		// For the case of existing words for this target taxonomic group,
		// we need to merge the original and new words, and then make them
		// unique
		//deque<Word> tmp;
		
		//merge(	w.begin(), w.end(), 
		//	local_w.begin(), local_w.end(),
		//	back_inserter(tmp) );
		
		//tmp.erase( unique( tmp.begin(), tmp.end() ), tmp.end() );
		
		//w.swap(tmp);
		unique_merge(w, local_w);
	}
	
	// Return the size of the number of target words we just added. Since the new words may be
	// a subset of the old words, return a size of a least sizeof(size_t) bytes to (a) account
	// for the size of the target index and (b) make sure this target is included in the list 
	// of valid targets (we use a size of 0 to indicate *invalid* targets).
	return (w.size() - init_word_count)*sizeof(Word) + sizeof(size_t);
}
