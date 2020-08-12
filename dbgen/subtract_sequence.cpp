#include "gottcha_db.h"

#include <queue>

#include <omp.h>

using namespace std;

void subtract_sequence(deque< pair< GroupInfo, deque<Word> > > &m_target_words, 
	const deque< pair< GroupInfo, deque<Word> > > &m_background_words,
	UpdateInfo &m_progress, bool m_update_progress, bool m_subtract_target)
{
	const size_t num_target_group = m_target_words.size();
	const size_t num_background_group = m_background_words.size();
	
	// DEBUG
	//size_t min_num_words = 0xFFFFFFFFFFFFFFFF;
	//size_t max_num_words = 0;
	
	#pragma omp parallel \
		firstprivate(m_update_progress)
	{
		const int num_thread = omp_get_num_threads();
		const int tid = omp_get_thread_num();
		
		// DEBUG
		//size_t num_thread_words = 0;
		
		// Only the first thread updates the user
		m_update_progress = m_update_progress && (tid == 0);
		
		priority_queue< pair<Word, size_t> > target_buffer;
		priority_queue< pair<Word, size_t> > background_buffer;

		vector< pair<deque<Word>::reverse_iterator, size_t> > target_scratch;
		vector<size_t> background_scratch;

		// Track the current word for each taxonomic group (starting at the back 
		// so that we don't need to change the behaviour of the priority_queue, which
		// returns the *largest* value).
		vector< deque<Word>::reverse_iterator > target_iter(num_target_group);
		vector< deque<Word>::const_reverse_iterator > background_iter(num_background_group);

		// For tracking progress
		size_t num_target_words = 0;

		for(size_t i = 0;i < num_target_group;++i){

			target_iter[i] = m_target_words[i].second.rbegin();
	
			while( (m_target_words[i].second.rend() != target_iter[i]) &&
				(kmer_hash(target_iter[i]->kmer & THREAD_HASH_BITS, num_thread) != tid) ){
				
				++target_iter[i];
			}
			
			if( m_target_words[i].second.rend() != target_iter[i] ){
				
				target_buffer.push( make_pair(*target_iter[i], i) );
				
				// DEBUG
				//++num_thread_words;
			}

			num_target_words += m_target_words[i].second.size();
		}

		// Scale the number of target words by the number of threads (assuming
		// an even mapping of words to threads)
		num_target_words = max( size_t(1), num_target_words/num_thread );
		
		for(size_t i = 0;i < num_background_group;++i){

			background_iter[i] = m_background_words[i].second.rbegin();
			
			while( (m_background_words[i].second.rend() != background_iter[i]) &&
				(kmer_hash(background_iter[i]->kmer & THREAD_HASH_BITS, num_thread) != tid) ){
				
				++background_iter[i];
			}
			
			if( m_background_words[i].second.rend() != background_iter[i] ){
				background_buffer.push( make_pair(*background_iter[i], i) );
			}
		}

		// Keep the user informed
		const size_t update_every = max( size_t(1), num_target_words/100 );
		size_t curr_num_target_words = 0;
		size_t next_update = update_every;

		while( !target_buffer.empty() && !background_buffer.empty() ){

			// The largest word in the priority queue
			const pair<Word, size_t> max_target_word = target_buffer.top();
			const pair<Word, size_t> max_background_word = background_buffer.top();

			if(max_target_word.first < max_background_word.first){

				// This word does *not* occur in any target.

				background_buffer.pop();

				// Are there additional words still waiting in the background taxa group?
				do{
					++background_iter[max_background_word.second];
				}
				while( (m_background_words[max_background_word.second].second.rend() != background_iter[max_background_word.second]) &&
					(kmer_hash(background_iter[max_background_word.second]->kmer & THREAD_HASH_BITS, num_thread) != tid) );

				if(m_background_words[max_background_word.second].second.rend() != 
					background_iter[max_background_word.second]){

					background_buffer.push( 
						make_pair(*background_iter[max_background_word.second], 
							max_background_word.second) );
				}

				continue;
			}

			//////////////////////////////////////////////////////////////////////////////
			// This word is found in *one* or more targets and *zero* or more backgrounds
			//////////////////////////////////////////////////////////////////////////////

			// Handle all of the *target* groups that contain max_target_word.
			while( !target_buffer.empty() && (target_buffer.top().first == max_target_word.first) ){

				const unsigned int index = target_buffer.top().second;

				target_buffer.pop();
				++curr_num_target_words;

				target_scratch.push_back( make_pair(target_iter[index], index) );
				
				do{
					++target_iter[index];
				}
				while( (m_target_words[index].second.rend() != target_iter[index]) &&
					(kmer_hash(target_iter[index]->kmer & THREAD_HASH_BITS, num_thread) != tid) );

				if(m_target_words[index].second.rend() != target_iter[index] ){
					
					target_buffer.push( make_pair(*target_iter[index], index) );
					
					// DEBUG
					//++num_thread_words;
				}
			}

			// Handle all of the *background* groups that contain *max_target_word.*
			while( !background_buffer.empty() && (background_buffer.top().first == max_target_word.first) ){

				const unsigned int index = background_buffer.top().second;

				background_buffer.pop();

				background_scratch.push_back(index);
				
				do{
					++background_iter[index];
				}
				while( (m_background_words[index].second.rend() != background_iter[index]) &&
					(kmer_hash(background_iter[index]->kmer & THREAD_HASH_BITS, num_thread) != tid) );
			
				if(m_background_words[index].second.rend() != background_iter[index] ){
					background_buffer.push( make_pair(*background_iter[index], index) );
				}
			}

			// Mask the target words
			const size_t num_target_scratch = target_scratch.size();

			for(size_t i = 0;i < num_target_scratch;++i){

				Word& w_a = *(target_scratch[i].first);
				const GroupInfo& group_a = m_target_words[target_scratch[i].second].first;

				// compare against background
				for(vector<size_t>::const_iterator j = background_scratch.begin();
					j != background_scratch.end();++j){

					const GroupInfo& group_b = m_background_words[*j].first;
					const TaxaSet mask = group_mask(group_a, group_b);

					w_a.taxa &= mask;
				}

				if(m_subtract_target){

					// compare against target
					for(size_t j = i + 1;j < num_target_scratch;++j){

						Word& w_b = *(target_scratch[j].first);
						const GroupInfo& group_b = m_target_words[target_scratch[j].second].first;

						const TaxaSet mask = group_mask(group_a, group_b);

						w_a.taxa &= mask;
						w_b.taxa &= mask;
					}
				}
			}

			// We must manually clean up scratch space after every iteration
			target_scratch.clear();
			background_scratch.clear();

			if( m_update_progress && (curr_num_target_words >= next_update) ){

				m_progress << "Subtracting words: " 
					<< 100.0*double(curr_num_target_words)/num_target_words 
					<< '%';
				m_progress.flush();

				next_update = (curr_num_target_words/update_every + 1)*update_every;
			}
		}

		// If we are subtracting targets from each other *and* there is remaining target
		// words, then we need to subtract them.
		while(m_subtract_target && !target_buffer.empty() ){

			// The largest word in the priority queue
			const pair<Word, size_t> max_target_word = target_buffer.top();

			// Handle all of the *target* groups that contain max_target_word.
			while( !target_buffer.empty() && (target_buffer.top().first == max_target_word.first) ){

				const unsigned int index = target_buffer.top().second;

				target_buffer.pop();
				++curr_num_target_words;

				target_scratch.push_back( make_pair(target_iter[index], index) );
				
				do{
					++target_iter[index];
				}
				while( (m_target_words[index].second.rend() != target_iter[index]) &&
					(kmer_hash(target_iter[index]->kmer & THREAD_HASH_BITS, num_thread) != tid) );
			
				if(m_target_words[index].second.rend() != target_iter[index] ){

					target_buffer.push( make_pair(*target_iter[index], index) );
					
					// DEBUG
					//++num_thread_words;
				}
			}

			// Mask the target words
			const size_t num_target_scratch = target_scratch.size();

			if(num_target_scratch > 1){

				for(size_t i = 0;i < num_target_scratch;++i){

					Word& w_a = *(target_scratch[i].first);
					const GroupInfo& group_a = m_target_words[target_scratch[i].second].first;

					// compare against target
					for(size_t j = i + 1;j < num_target_scratch;++j){

						Word& w_b = *(target_scratch[j].first);
						const GroupInfo& group_b = m_target_words[target_scratch[j].second].first;

						const TaxaSet mask = group_mask(group_a, group_b);

						w_a.taxa &= mask;
						w_b.taxa &= mask;
					}
				}
			}

			// We must manually clean up scratch space after every iteration
			target_scratch.clear();

			if( m_update_progress && (curr_num_target_words >= next_update) ){

				m_progress << "Subtracting words: " 
					<< 100.0*double(curr_num_target_words)/num_target_words 
					<< '%';
				m_progress.flush();

				next_update = (curr_num_target_words/update_every + 1)*update_every;
			}
		}
		
		//#pragma omp critical
		//{
		//	min_num_words = min(min_num_words, num_thread_words);
		//	max_num_words = max(max_num_words, num_thread_words);
		//}
	}
	
	// DEBUG
	//cerr << "\nSubtraction balance [" << mpi_rank << "] = " 
	//	<< double(min_num_words)/max_num_words << endl;
}
