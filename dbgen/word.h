#ifndef __WORD
#define __WORD

// Indicate the set of taxa levels that a word is unique to
// (one bit per taxa level, so if the list of taxa levels in
// options.h grows to be larger than eight, we will need to
// increase the size of TaxaSet to an unsigned short).
typedef unsigned char TaxaSet;
#define	NO_TAXA			0x00

typedef size_t BaseWord;

// Using two bits per base, we can fit a 58 bit (28 bases) kmer and
// 8 bits of taxa set into a 64-bit size_t.
struct Word
{
	BaseWord kmer : 56; // A packed kmer
	TaxaSet taxa : 8;  // Bits indicating the set of taxa this kmer is unqiue to
		
	//BaseWord kmer; // A packed kmer
	//TaxaSet taxa;  // Bits indicating the set of taxa this kmer is unqiue to
	
	Word()
	{
		kmer = 0;
		taxa = 0;
	};
	
	// A single argument to the constructor initializes the kmer
	Word(const BaseWord &m_kmer) :
		kmer(m_kmer)
	{
		taxa = 0;
	};
	
	Word(const BaseWord &m_kmer, const TaxaSet &m_taxa) :
		kmer(m_kmer), taxa(m_taxa)
	{
	};
		
	// Comparisons use the kmer portion
	inline bool operator<(const Word &m_rhs) const
	{
		return (kmer < m_rhs.kmer);
	}
	
	inline bool operator>(const Word &m_rhs) const
	{
		return (kmer > m_rhs.kmer);
	}
	
	inline bool operator==(const Word &m_rhs) const
	{
		return (kmer == m_rhs.kmer);
	}
	
	inline bool operator!=(const Word &m_rhs) const
	{
		return (kmer != m_rhs.kmer);
	}
};

// Add an index to a word
struct IndexedWord
{
	BaseWord kmer;
	size_t index;
	
	IndexedWord(const BaseWord &m_base_word, const size_t &m_index) :
		kmer(m_base_word), index(m_index)
	{
	};
	
	// Comparisons use the kmer as the primary sorting key and the
	// index as the secondary sorting key
	inline bool operator<(const IndexedWord &m_rhs) const
	{
		return (kmer == m_rhs.kmer) ?
			(index < m_rhs.index) :
			(kmer < m_rhs.kmer);
	}
};

inline int kmer_hash(BaseWord m_kmer, const int &m_num_task)
{
	// From http://burtleburtle.net/bob/hash/integer.html
	// This website has a number of pretty good hash functions.
	// Even though these functions expect 32 bit input, they still seem
	// to work well for the (up to) 64 bit inputs used below (I've currently
	// only tested with 48 bit Words).
	// Hash #1
	//w = (w^0xdeadbeef) + (w<<4);
	//w = w ^ (w>>10);
	//w = w + (w<<7);
	//w = w ^ (w>>13);
	//return w % m_num_tasks;

	// Hash #2
	//w = w ^ (w>>4);
    	//w = (w^0xdeadbeef) + (w<<5);
    	//w = w ^ (w>>11);
    	//return w % m_num_tasks;

	// Hash #3
    	m_kmer -= (m_kmer << 6);
    	m_kmer ^= (m_kmer >> 17);
    	m_kmer -= (m_kmer << 9);
    	m_kmer ^= (m_kmer << 4);
    	m_kmer -= (m_kmer << 3);
    	m_kmer ^= (m_kmer << 10);
    	m_kmer ^= (m_kmer >> 15);
	
    	return m_kmer % m_num_task;
}

#define MAX_WORD_SIZE	28
#endif // __WORD
