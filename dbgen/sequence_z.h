#ifndef __SEQUENCE_Z
#define __SEQUENCE_Z

#include <string>
#include <deque>
#include "mpi_util.h"

class SequenceZ {
	public:
		enum {A, G, C, T, INVALID_BASE};
    
	private:
    
		//typedef size_t Element;
		typedef unsigned int Element;
    
    		// Two bits per base, 4 bases per byte
		#define     BASE_PER_ELEM   ( sizeof(SequenceZ::Element)*4 )

		struct Fragment : public std::deque<Element> {
            		size_t start;
            		size_t len; // The number of bases packed into the deque
            
            		Fragment(const size_t &m_start)
			{
				start = m_start;
				len = 0;
			}
			
			inline void pack(unsigned char m_base, unsigned char &m_offset)
			{
				if(m_offset == BASE_PER_ELEM){
					
					m_offset = 0;
					push_back( Element() );
				}
				
				back() |= Element(m_base) << 2*m_offset;
				
				++len;
				++m_offset;
			};
        	};
    
        	std::deque<Fragment> data;

        	// The offset is only used for pushing bases onto the growing sequence
        	unsigned char offset;

        	size_t min_fragment_size;
		
	public:

		SequenceZ(const size_t &m_min_fragment_size = 0)
		{
			offset = BASE_PER_ELEM;
			min_fragment_size = m_min_fragment_size;
			
			// Create an empty fragment that starts at 0
			data.push_back( Fragment(0) );
		};
	
        	class const_iterator{
        	    private:
                	std::deque<Fragment>::const_iterator iter;
                	std::deque<Fragment>::const_iterator end_iter;
                	Fragment::const_iterator frag_iter;
			size_t num_frag_base;
                	unsigned char base_offset;
			bool is_new_frag;
			
        	    public:

			const_iterator(
			    const std::deque<Fragment>::const_iterator &m_begin_iter,
			    const std::deque<Fragment>::const_iterator &m_end_iter) :
        			iter(m_begin_iter), end_iter(m_end_iter)
			{
				num_frag_base = 0;
				base_offset = 0;
				is_new_frag = true;

				if(m_begin_iter != m_end_iter){
					frag_iter = iter->begin();
				}
			};

			const_iterator(const const_iterator &m_rhs) :
				iter(m_rhs.iter),
				end_iter(m_rhs.end_iter),
				frag_iter(m_rhs.frag_iter),
				num_frag_base(m_rhs.num_frag_base),
				base_offset(m_rhs.base_offset),
				is_new_frag(m_rhs.is_new_frag)
			{
			};
			
			inline unsigned char binary_base() const
			{
				return (*frag_iter >> 2*base_offset) & 3;
			};
			
			inline char operator*() const
			{
				switch( binary_base() ){
					case SequenceZ::A:
						return 'A';
					case SequenceZ::T:
						return 'T';
					case SequenceZ::G:
						return 'G';
					case SequenceZ::C:
						return 'C';
					default:
						throw ":SequenceZ::const_iterator::operator*: Unknown base!";
				};

				// We should never get here!
				return '?';
			};

                	inline bool operator!=(const const_iterator &m_rhs) const
                	{
				// Note that we only test "frag_iter != m_rhs.frag_iter"
				// after checking that "iter != end_iter" is true (since
				// m_rhs.frag_iter is not initialized for end_iter.
				return (iter != m_rhs.iter) ||
					( (iter != end_iter) && (frag_iter != m_rhs.frag_iter) ) ||
					(base_offset != m_rhs.base_offset);
                	};

			inline bool operator==(const const_iterator &m_rhs) const
                	{
                	    return !(*this != m_rhs);
                	};
			
                	// Prefix (only increment allowed -- not postfix operator provided!)
                	inline const_iterator& operator++()
			{
				++num_frag_base;
				++base_offset;
				is_new_frag = false;

				if(base_offset == BASE_PER_ELEM){

					++frag_iter;
					base_offset = 0;
				}

				// Sequences can have empty terminal fragments with len = 0
				// so use "<=" in the following test:
				if(iter->len <= num_frag_base){
				
					++iter;
					num_frag_base = 0;
					base_offset = 0;
					is_new_frag = true;

					if(iter != end_iter){
						frag_iter = iter->begin();
					}
				}
				
				return *this;
			};

			inline const_iterator operator+(size_t m_i) const
			{
				const_iterator ret(*this);
				
				//// The naive way ...
				//for(size_t i = 0;i < m_i;++i){
				//	std::cerr << "i = " << i << std::endl;
				//
				//	++ret;
				//}
				//return ret;
				
				// Not allowed to increment past the
				// end of the sequence
				if(ret.iter == ret.end_iter){
					
					if(m_i == 0){
						return ret;
					}
					
					throw __FILE__ ":SequenceZ::const_iterator+: Can't increment past the end (1)";
				}
				
				// Move to the end of the current fragment if we can
				while( m_i > (ret.iter->len - ret.num_frag_base) ){
					
					m_i -= (ret.iter->len - ret.num_frag_base);
					
					++ret.iter;
					ret.num_frag_base = 0;
					ret.base_offset = 0;
					ret.is_new_frag = true;
					
					if(ret.iter == ret.end_iter){
						throw __FILE__ ":SequenceZ::const_iterator+: Can't increment past the end (2)";
					}
					
					ret.frag_iter = ret.iter->begin();
				}
				
				// Move to the next element boundary if we can
				while( (m_i > 0) && (ret.base_offset != 0) ){
					
					--m_i;
					++ret.num_frag_base;
					++ret.base_offset;
					ret.is_new_frag = false;
					
					if(ret.base_offset == BASE_PER_ELEM){
						
						++ret.frag_iter;
						ret.base_offset = 0;
					}
				}
				
				const size_t delta_elem = m_i/BASE_PER_ELEM;
				const size_t delta_base = delta_elem*BASE_PER_ELEM;

				m_i -= delta_base;
				ret.num_frag_base += delta_base;
				ret.frag_iter += delta_elem;
				
				// The remainder
				while(m_i > 0){
					
					--m_i;
					++ret.num_frag_base;
					++ret.base_offset;
					ret.is_new_frag = false;
					
					if(ret.base_offset == BASE_PER_ELEM){
						
						++ret.frag_iter;
						ret.base_offset = 0;
					}
				}
				
				// Sequences can have empty terminal fragments with len = 0
				// so use "<=" in the following test:
				if(ret.iter->len <= ret.num_frag_base){
				
					++ret.iter;
					ret.num_frag_base = 0;
					ret.base_offset = 0;
					ret.is_new_frag = true;

					if(ret.iter != ret.end_iter){
						ret.frag_iter = ret.iter->begin();
					}
				}
				
				return ret;
			};
			
                	inline size_t index () const
			{
				return iter->start + num_frag_base;
			};
			
			inline bool new_fragment() const
			{
				return is_new_frag;
			};
        	};

		inline void push_back(char m_base){
			
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
		};
	
        	inline const_iterator begin() const
        	{
        	    // data should never be empty, but it can contain an empty fragment.
        	    if(data.front().len == 0){

                	    // Skip the empty fragment at the begining of data
                	    return const_iterator( data.begin() + 1, data.end() );
        	    }

        	    return const_iterator( data.begin(), data.end() );
        	};

        	inline const_iterator end() const
        	{
        	    return const_iterator( data.end(), data.end() );
        	};

		inline void clear()
		{
			offset = BASE_PER_ELEM;
			
			data.clear();
			
			// Create an empty fragment that starts at 0
			data.push_back( Fragment(0) );
		};

        	inline size_t size() const
        	{
        	    size_t ret = 0;

        	    for(std::deque<Fragment>::const_iterator i = data.begin();
                	i != data.end();++i){
                	ret += i->len;
        	    }

        	    return ret;
        	};

		// Circularize a sequence by adding m_len - 1bases from the
		// front to the back
		void circularize(const size_t &m_len)
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
		};

        	template<class T> friend
        	    size_t mpi_size(const T &m_obj);
        	template<class T> friend
        	    unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
        	template<class T> friend
        	   unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const SequenceZ &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceZ &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceZ &m_obj);

#endif // __SEQUENCE_Z
