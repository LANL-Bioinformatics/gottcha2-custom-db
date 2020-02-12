#include "mpi_util.h"
#include "gottcha_db.h"
#include "sequence_z.h"
#include <string.h>

using namespace std;

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for std::string
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<string>(const string &m_str)
{
	return sizeof(size_t) + m_str.size();
}

template<>
unsigned char* mpi_pack<string>(unsigned char* m_ptr, const string &m_str)
{
	size_t len = m_str.size();
	
	memcpy( m_ptr, &len, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	memcpy(m_ptr, m_str.c_str(), len);
	m_ptr += len;
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<string>(unsigned char* m_ptr, string &m_str)
{
	size_t len;
	
	memcpy( &len, m_ptr, sizeof(size_t) );
	m_ptr += sizeof(size_t);
	
	m_str.assign( (char*)m_ptr, len );
	m_ptr += len;
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for __uint128_t
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size<__uint128_t>(const __uint128_t &m_obj)
{
	return sizeof(__uint128_t);
}

template<>
unsigned char* mpi_pack<__uint128_t>(unsigned char* m_ptr, const __uint128_t &m_obj)
{
	memcpy( m_ptr, &m_obj, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack<__uint128_t>(unsigned char* m_ptr, __uint128_t &m_obj)
{
	memcpy( &m_obj, m_ptr, sizeof(__uint128_t) );
	m_ptr += sizeof(__uint128_t);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for GroupInfo
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size(const GroupInfo &m_obj)
{
	return NUM_TAXA_LEVEL*sizeof(unsigned int);
}

template<>
unsigned char* mpi_pack(unsigned char* m_ptr, const GroupInfo &m_obj)
{
		
	memcpy( m_ptr, m_obj.id, NUM_TAXA_LEVEL*sizeof(unsigned int) );
	m_ptr += NUM_TAXA_LEVEL*sizeof(unsigned int);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack(unsigned char* m_ptr, GroupInfo &m_obj)
{
	memcpy( m_obj.id, m_ptr, NUM_TAXA_LEVEL*sizeof(unsigned int) );
	m_ptr += NUM_TAXA_LEVEL*sizeof(unsigned int);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SequenceFileInfo
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size(const SequenceFileInfo &m_obj)
{
	return mpi_size(m_obj.taxa_id) + mpi_size(m_obj.topo) + 
		mpi_size(m_obj.group_id);
}

template<>
unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceFileInfo &m_obj)
{
		
	m_ptr = mpi_pack(m_ptr, m_obj.taxa_id);
	m_ptr = mpi_pack(m_ptr, m_obj.topo);
	m_ptr = mpi_pack(m_ptr, m_obj.group_id);
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceFileInfo &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.taxa_id);
	m_ptr = mpi_unpack(m_ptr, m_obj.topo);
	m_ptr = mpi_unpack(m_ptr, m_obj.group_id);
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SequenceZ
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size(const SequenceZ &m_obj)
{
	size_t ret = mpi_size(m_obj.offset) + 
		mpi_size(m_obj.min_fragment_size) +
		sizeof(size_t); // Number of fragments
	
	for(std::deque<SequenceZ::Fragment>::const_iterator i = m_obj.data.begin();
		i != m_obj.data.end();++i){
		
		// Number of elements in fragments, start and len
		ret += 3*sizeof(size_t) +
			i->size()*sizeof(SequenceZ::Element);
	}
	
	return ret;
}

template<>
unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceZ &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.offset);
	m_ptr = mpi_pack(m_ptr, m_obj.min_fragment_size);
	
	const size_t num_frag = m_obj.data.size();
	
	m_ptr = mpi_pack(m_ptr, num_frag);
	
	for(std::deque<SequenceZ::Fragment>::const_iterator i = m_obj.data.begin();
		i != m_obj.data.end();++i){
		
		const size_t num_elem = i->size();
		
		m_ptr = mpi_pack(m_ptr, num_elem);
		m_ptr = mpi_pack(m_ptr, i->start);
		m_ptr = mpi_pack(m_ptr, i->len);
		
		for(size_t j = 0;j < num_elem;++j){
			m_ptr = mpi_pack(m_ptr, (*i)[j]);
		}
	}
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceZ &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.offset);
	m_ptr = mpi_unpack(m_ptr, m_obj.min_fragment_size);
	
	size_t num_frag = 0;
	
	m_ptr = mpi_unpack(m_ptr, num_frag);
	
	m_obj.data.clear();
	
	for(size_t i = 0;i < num_frag;++i){
		
		m_obj.data.push_back( SequenceZ::Fragment(0) );
		
		SequenceZ::Fragment &ref = m_obj.data.back();
		
		size_t num_elem = 0;
		
		m_ptr = mpi_unpack(m_ptr, num_elem);
		m_ptr = mpi_unpack(m_ptr, ref.start);
		m_ptr = mpi_unpack(m_ptr, ref.len);
		
		for(size_t j = 0;j < num_elem;++j){
			
			SequenceZ::Element local;
			
			m_ptr = mpi_unpack( m_ptr, local);
			
			ref.push_back(local);
		}
	}
	
	return m_ptr;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for SequenceFile
/////////////////////////////////////////////////////////////////////////////////////////
template<>
size_t mpi_size(const SequenceFile &m_obj)
{
	size_t ret = mpi_size(m_obj.info) + mpi_size(m_obj.name) + sizeof(size_t);
	
	for(deque< pair<string, SequenceZ*> >::const_iterator i = m_obj.seq.begin();
		i != m_obj.seq.end();++i){
		
		ret += mpi_size(i->first) + mpi_size( *(i->second) );
	}
	
	return ret;
}

template<>
unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceFile &m_obj)
{
	m_ptr = mpi_pack(m_ptr, m_obj.info);
	m_ptr = mpi_pack(m_ptr, m_obj.name);
	
	const size_t num_seq = m_obj.seq.size();
	
	m_ptr = mpi_pack(m_ptr, num_seq);
	
	for(deque< pair<string, SequenceZ*> >::const_iterator i = m_obj.seq.begin();
		i != m_obj.seq.end();++i){
		
		m_ptr = mpi_pack(m_ptr, i->first);
		m_ptr = mpi_pack(m_ptr, *(i->second) );
	}
	
	return m_ptr;
}

template<>
unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceFile &m_obj)
{
	m_obj.clear();
	
	m_ptr = mpi_unpack(m_ptr, m_obj.info);
	m_ptr = mpi_unpack(m_ptr, m_obj.name);

	size_t num_seq;
	
	m_ptr = mpi_unpack(m_ptr, num_seq);
	
	m_obj.seq.resize(num_seq);
	
	for(size_t i = 0;i < num_seq;++i){
		
		m_ptr = mpi_unpack(m_ptr, m_obj.seq[i].first);
		
		m_obj.seq[i].second = new SequenceZ;
		
		if(m_obj.seq[i].second == NULL){
			throw __FILE__ ":mpi_unpack<SequenceFile>: Unable to allocate sequence";
		}
		
		m_ptr = mpi_unpack(m_ptr, *( m_obj.seq[i].second ) );
	}
	
	return m_ptr;
}

/////////////////////////////////////////////////////////////////////////////////////////
// Specialization for Options
/////////////////////////////////////////////////////////////////////////////////////////

template<> 
size_t mpi_size(const Options &m_obj)
{
	return  mpi_size(m_obj.word_size) + 
		mpi_size(m_obj.min_output_frag) + 
		mpi_size(m_obj.compress_output) + 
		mpi_size(m_obj.quit) + 
		mpi_size(m_obj.default_topology) +
		mpi_size(m_obj.output_root_dir) +
		mpi_size(m_obj.taxa_level_output_prefix) +
		
		mpi_size(m_obj.taxa_level_mapping_file) +
		mpi_size(m_obj.groups_to_squash) +
		mpi_size(m_obj.log_file) +
		mpi_size(m_obj.max_worker_RAM);
};

template<> 
unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj)
{
	m_ptr = mpi_unpack(m_ptr, m_obj.word_size);
	m_ptr = mpi_unpack(m_ptr, m_obj.min_output_frag);
	m_ptr = mpi_unpack(m_ptr, m_obj.compress_output);
	m_ptr = mpi_unpack(m_ptr, m_obj.quit);
	m_ptr = mpi_unpack(m_ptr, m_obj.default_topology);
	m_ptr = mpi_unpack(m_ptr, m_obj.output_root_dir);
	m_ptr = mpi_unpack(m_ptr, m_obj.taxa_level_output_prefix);
	m_ptr = mpi_unpack(m_ptr, m_obj.taxa_level_mapping_file);
	m_ptr = mpi_unpack(m_ptr, m_obj.groups_to_squash);
	m_ptr = mpi_unpack(m_ptr, m_obj.log_file);
	m_ptr = mpi_unpack(m_ptr, m_obj.max_worker_RAM);
	
	return m_ptr;
}

template<> 
unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj)
{

	m_ptr = mpi_pack(m_ptr, m_obj.word_size);
	m_ptr = mpi_pack(m_ptr, m_obj.min_output_frag);
	m_ptr = mpi_pack(m_ptr, m_obj.compress_output);
	m_ptr = mpi_pack(m_ptr, m_obj.quit);
	m_ptr = mpi_pack(m_ptr, m_obj.default_topology);
	m_ptr = mpi_pack(m_ptr, m_obj.output_root_dir);
	m_ptr = mpi_pack(m_ptr, m_obj.taxa_level_output_prefix);
	m_ptr = mpi_pack(m_ptr, m_obj.taxa_level_mapping_file);
	m_ptr = mpi_pack(m_ptr, m_obj.groups_to_squash);
	m_ptr = mpi_pack(m_ptr, m_obj.log_file);
	m_ptr = mpi_pack(m_ptr, m_obj.max_worker_RAM);

	return m_ptr;
}
