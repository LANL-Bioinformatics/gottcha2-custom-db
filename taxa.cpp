#include "gottcha_db.h"
#include "deque_set.h"

using namespace std;

TaxaSet get_taxa_mask(const deque<unsigned int> &m_groups_to_squash, const SequenceFileInfo &m_info)
{
	TaxaSet ret = NO_TAXA;

	for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

		if(m_info.group_id[level] == UNDEFINED_GROUP){

			// No computation is specified for this taxa level
			continue;
		}

		// There *may* be computation required for this taxa level. Set the level
		// bit now (but it may be cleared if this file is squashed at this, or a 
		// higher, taxonomic level).
		ret |= 1 << level;

		if(set_contains(m_groups_to_squash, m_info.group_id[level]) == true){

			// Clear all of the bits at the current and lower taxonomic
			// levels.
			ret = NO_TAXA;
		}		
	}
	
	return ret;
}

TaxaSet group_mask(const GroupInfo &m_background_info, const GroupInfo &m_target_info)
{
	TaxaSet ret = 0; // Start with no bits set
	
	for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){
		
		// Set a bit at each level where the taxonomic groups are the *same*
		ret |= ( (m_target_info[level] == m_background_info[level]) << level);
	}
	
	// Since the returned value will be used as a bit mask (using the "&" operator),
	// bits that are "1" (which correspond to the *same* taxonomic group at a given 
	// taxonomic level) will *not* change the membership state of a given word at
	// the specified taxonomic level. Different taxonomic groups at a given level
	// will yield a "0" in the mask and will force the word state to "0" at the 
	// same level.
	return ret;
}
