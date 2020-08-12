#include "taxa_id.h"

unsigned int num_digits(uint64_t m_number);

TaxaId str_to_taxa_id(const std::string &m_tid, bool m_ignore_strain /*= false*/)
{
	uint64_t taxa_id = 0;
	uint64_t strain_id = 0;
	
	// The first digits we read are the NCBI taxa id
	uint64_t *curr = &taxa_id;
	
	for(std::string::const_iterator i = m_tid.begin();i != m_tid.end();++i){
		
		switch(*i){
			case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
			case '8': case '9':
				*curr = (*curr)*10 + (*i - '0');
				break;
			case '.':
				if(curr == &strain_id){
					throw __FILE__ ":str_to_taxa_id: Multiple '.' found!";
				}
				
				// Start reading the strain id
				curr = &strain_id;
				
				break;
			default:
				throw __FILE__ ":str_to_taxa_id: Invalid symbol";
		};
	}

	return m_ignore_strain ? 
		(taxa_id << 32) :
		(taxa_id << 32) | (strain_id);
}

std::string taxa_id_to_str(const TaxaId &m_tid, bool m_show_strain /*= true*/)
{
	std::string ret;
	uint64_t taxa_id = (m_tid >> 32) & 0xFFFFFFFF;
	uint64_t strain_id = m_tid & 0xFFFFFFFF;

	unsigned int index = num_digits(taxa_id) + num_digits(strain_id) + (m_show_strain ? 1 : 0);
	
	ret.resize(index);
	
	if(m_show_strain){
		
		do{
		
			ret[--index] = (strain_id % 10) + '0';
			strain_id /= 10;
		}
		while(strain_id > 0);
		
		ret[--index] = '.';
	}
	
	do{
	
		ret[--index] = (taxa_id % 10) + '0';
		taxa_id /= 10;
	}
	while(taxa_id > 0);
	
	return ret;
}

unsigned int num_digits(uint64_t m_number)
{
	unsigned int ret = (m_number == 0) ? 1 : 0;
	
	while(m_number > 0){
		
		++ret;
		m_number /= 10;
	}
	
	return ret;
}
