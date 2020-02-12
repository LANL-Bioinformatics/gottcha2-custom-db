#include "sequence_z.h"
#include "mpi_util.h"

SequenceZ::SequenceZ(const size_t &m_min_fragment_size)
{
	offset = BASE_PER_ELEM;
	min_fragment_size = m_min_fragment_size;

	// Create an empty fragment that starts at 0
	data.push_back( Fragment(0) );
};

void SequenceZ::clear()
{
	offset = BASE_PER_ELEM;

	data.clear();

	// Create an empty fragment that starts at 0
	data.push_back( Fragment(0) );
}

SequenceZ::Fragment::Fragment(const size_t &m_start)
{
	start = m_start;
	len = 0;
}

void SequenceZ::Fragment::pack(unsigned char m_base, unsigned char &m_offset)
{
	if(m_offset == BASE_PER_ELEM){
        
		m_offset = 0;
		push_back( Element() );
	}

	back() |= Element(m_base) << 2*m_offset;

	++len;
	++m_offset;
}

void SequenceZ::push_back(char m_base)
{
	if( data.empty() ){
		throw __FILE__ ":SequenceZ::push_back: Empty fragment buffer";
	}

	Fragment &f = data.back();

	switch(m_base){
		case 'A': case 'a':
			f.pack(SequenceZ::A, offset);
			break;
		case 'T': case 't':
			f.pack(SequenceZ::T, offset);
			break;
		case 'G': case 'g':
			f.pack(SequenceZ::G, offset);
			break;
		case 'C': case 'c':
			f.pack(SequenceZ::C, offset);
			break;
		case ' ': case '\t': case '\r': case '\n':
			// Ignore white space
			break;
		default:
			
			// A non-ATGC base, like N
			if(f.len < min_fragment_size){

				f.clear();
				f.start = f.len + f.start + 1;
				f.len = 0;
			}
			else{
				data.push_back( Fragment(f.len + f.start + 1) );
			}

			offset = BASE_PER_ELEM;

			break;
	};
}

// Add m_len - 1 bases from the front of the sequence to the back
void SequenceZ::circularize(const size_t &m_len)
{
	if(m_len == 0){
		return;
	}
	
	std::deque<char> buffer;
	size_t index = 0;
	
	for(const_iterator i = begin();i != end();++i){
		
		buffer.push_back(*i);
		
		++index;
		
		if( (index == (m_len - 1) ) ||
		    ( (index > 1) && i.new_fragment() ) ){
			break;
		}
	}
	
	for(std::deque<char>::const_iterator i = buffer.begin();
		i != buffer.end();++i){
		
		push_back(*i);
	}
}
