#ifndef __GOTTCHA_DB
#define __GOTTCHA_DB

#include <string>
#include <deque>
#include <vector>
#include <list>
#include <map>
#include <mpi.h>

#include "word.h"
#include "options.h"
#include "update.h"
#include "sequence_z.h"
#include "mpi_util.h"

#define	VERSION	"2.85"

// A worker RAM value of 0 is a special value that directs each
// compute node to measure the amount of available RAM prior to
// initiating the signature calculations and use at most 
// RAM_FRACTION of this amount.
#define		RAM_FRACTION		0.8

#define		UNDEFINED_GROUP		0

#define	PATH_SEPARATOR			'/'

#define	SEC_PER_MIN			60.0
#define	KB				1024
#define	MB				(1024*KB)
#define	GB				(1024*MB)

// For each base, there are two possible "sharing" states:
enum {
	YES = 0, // Shared
	NO,      // Not shared
};

// Set the NO bit for each of the eight possible taxonomic levels
#define	ALL_NO	0xFF

#define	RANK_HASH_BITS			0xAAAAAAAAAAAAAAAA
#define	THREAD_HASH_BITS		0x5555555555555555

// Sequence file types
typedef enum {
	FASTA, 
	GENBANK, 
	UNKNOWN_FILETYPE} FileType;


// MPI messages
enum {
	GOTTCHA_QUIT = 1000,
	GOTTCHA_FILE_SIZE,
	GOTTCHA_FILE_LOAD,
	GOTTCHA_WORD_SIZE,
	GOTTCHA_TARGET_FILE_REQUEST,
	GOTTCHA_BACKGROUND_FILE_REQUEST,
	GOTTCHA_PROFILE_LOAD,
	GOTTCHA_PROFILE_COMPUTE,
	GOTTCHA_PROFILE_TARGET_UPDATE,
	GOTTCHA_HOSTNAME
};

//////////////////////////////////////////////////////////////////////////////////

class GroupInfo
{
	private:
		unsigned int id[NUM_TAXA_LEVEL];
	
	public:
		unsigned int& operator[](const size_t &m_index)
		{
			return id[m_index];
		};
		
		unsigned int operator[](const size_t &m_index) const
		{
			return id[m_index];
		};
		
		GroupInfo()
		{
			memset( id, UNDEFINED_GROUP, NUM_TAXA_LEVEL*sizeof(unsigned int) );
		};
		
		inline bool operator<(const GroupInfo &m_rhs) const
		{
			for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

				if(id[level] != m_rhs.id[level]){
					return (id[level] < m_rhs.id[level]);
				}
			}
			
			return false;
		}
		
		inline bool operator==(const GroupInfo &m_rhs) const
		{
			for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

				if(id[level] != m_rhs.id[level]){

					// Return as soon as we find a *difference* in the group id
					return false;
				}
			}

			return true;
		};

		inline bool operator!=(const GroupInfo &m_rhs) const
		{
			for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){

				if(id[level] != m_rhs.id[level]){

					// Return as soon as we find a *difference* in the group id
					return true;
				}
			}

			return false;
		};
		
		template<class T> friend 
			size_t mpi_size(const T &m_obj);
		template<class T> friend 
			unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
		template<class T> friend 
			unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const GroupInfo &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const GroupInfo &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, GroupInfo &m_obj);

typedef std::string TaxaId;

struct SequenceFileInfo
{
	TaxaId taxa_id;
	Topology topo;
	GroupInfo group_id;
	
	SequenceFileInfo()
	{
		clear();
	};
	
	inline void clear()
	{
		topo = UNKNOWN_TOPOLOGY;
	}
	
	inline bool operator<(const SequenceFileInfo &m_rhs) const
	{
		// Sort by group id
		return (group_id < m_rhs.group_id);
	};
	
	template<class T> friend 
		size_t mpi_size(const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const SequenceFileInfo &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceFileInfo &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceFileInfo &m_obj);

namespace std 
{
	// Allow the use of GroupInfo as a hash key
	template <> struct hash<GroupInfo>
	{
		size_t operator()(const GroupInfo &m_info) const
		{
			size_t ret = 0;
			
			for(unsigned int level = 0;level < NUM_TAXA_LEVEL;++level){
				ret |= size_t(m_info[level]) << (level%2)*32;
			}
			
			return ret;
		};
	};
}

struct SequenceFile
{
	std::string name;
	SequenceFileInfo info;
	std::deque< std::pair<std::string, SequenceZ*> > seq;
	
	void clear()
	{
		for(std::deque< std::pair<std::string, SequenceZ*> >::iterator i = seq.begin();
			i != seq.end();++i){
			
			if(i->second != NULL){
				delete i->second;
			}
		}
		
		seq.clear();
		name.clear();
		info.clear();
	};
	
	// Return the sequence length
	inline size_t size() const
	{
		size_t ret = 0;
		
		for(std::deque< std::pair<std::string, SequenceZ*> >::const_iterator i = seq.begin();
			i != seq.end();++i){
			
			if(i->second != NULL){
				ret += i->second->size();
			}
		}
		
		return ret;
	};
	
	template<class T> friend 
		size_t mpi_size(const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
	template<class T> friend 
		unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const SequenceFile &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const SequenceFile &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, SequenceFile &m_obj);

//////////////////////////////////////////////////////////////////////////////////

// parse_mapping_file.cpp
void parse_mapping_files(
	std::vector<std::pair<SequenceFileInfo, std::string> > &m_file_info, 
	std::deque<unsigned int> &m_groups_to_squash, const Options &m_opt, 
	UpdateInfo &m_info);

// In io.cpp
void load_sequence(std::deque< std::pair<std::string, SequenceZ*> > &m_sequence, 
	const std::string &m_filename, const SequenceFileInfo &m_file_info,
	const size_t& m_word_size);
std::string remove_extension(const std::string &m_name, 
	const std::string &m_ext);
std::deque<std::string> split(const std::string &m_path, 
	const char &m_delim);	
FileType GetFileType(const std::string &m_path);
std::string load_sequence_database(
	std::unordered_map<size_t, std::pair<unsigned char*, size_t> > &m_db,
	std::vector<int> &m_sequence_location,
	std::deque<unsigned int> &m_groups_to_squash, 
	UpdateInfo &m_progress, const Options &m_opt);
void request_file(const size_t &m_file_index, const int &m_loc, 
	SequenceFile &m_background, 
	const std::unordered_map<size_t, std::pair<unsigned char*, size_t> > &m_db);
std::string format_output_filename(const std::string &m_input_filename, 
	const std::string &m_prefix, const std::string &m_root_dir, 
	const std::string &m_input_lcp, const bool &m_compress);
void create_directory(const std::string &m_path);
void create_complete_path(const std::string &m_filename);
std::string truncate_filename(const std::string &m_name, 
	const size_t &m_max_len);
size_t per_rank_database_size(const std::unordered_map<size_t, 
	std::pair<unsigned char*, size_t> > &m_db);

// In fasta.cpp
void parse_fasta(const std::string &m_filename, 
	std::deque< std::pair<std::string, SequenceZ*> > &m_sequence, 
	const SequenceFileInfo &m_file_info, const size_t& m_word_size);
std::pair<std::string, std::string> split_defline(const std::string &m_defline);
Topology get_topology(const SequenceFileInfo &m_file_info);
	
// In genbank.cpp
void parse_genbank(const std::string &m_filename, 
	std::deque< std::pair<std::string, SequenceZ*> > &m_sequence, 
	const SequenceFileInfo &m_file_info, const size_t& m_word_size);

// In meory_util.cpp
size_t available_RAM();
bool room_for_more(const bool &m_has_room);

// In write_targets.cpp
void write_targets(const std::vector<bool> &m_targets, 
	const std::deque< std::pair< GroupInfo, std::deque<Word> > > &m_target_words,
	const std::vector<int> &m_sequence_location,
	const std::unordered_map<size_t, std::pair<unsigned char*, size_t> > &m_db,
	UpdateInfo &m_progress, const std::string &m_input_lcp, 
	const Options &m_opt);

// In digest.cpp
void digest_words(std::deque<Word> &m_words, const SequenceZ &m_seq, 
	const unsigned int &m_word_size, const TaxaSet &m_init_taxa = NO_TAXA);
void digest_indexed_words(std::deque<IndexedWord> &m_words, const SequenceZ &m_seq, 
	const unsigned int &m_word_size, const size_t &m_index);	
size_t request_file_words(const size_t &m_target_index, const int &m_target_rank,
	std::string &m_target_name,
	std::deque< std::pair< GroupInfo, std::deque<Word> > > &m_target_words,
	const std::deque<unsigned int> &m_groups_to_squash, 
	const std::unordered_map<size_t, std::pair<unsigned char*, size_t> > &m_db,
	const size_t &m_word_size);

// In taxa.cpp
TaxaSet get_taxa_mask(const std::deque<unsigned int> &m_groups_to_squash, 
	const SequenceFileInfo &m_info);
TaxaSet group_mask(const GroupInfo &m_background_info, const GroupInfo &m_target_info);

// In subtract_sequence.cpp
void subtract_sequence(std::deque< std::pair< GroupInfo, std::deque<Word> > > &m_target_words, 
	const std::deque< std::pair< GroupInfo, std::deque<Word> > > &m_background_words,
	UpdateInfo &m_progress, bool m_update_progress, bool m_subtract_target);

#endif // __GOTTCHA_DB
