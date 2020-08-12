#ifndef __GOTTCHA_OPTIONS
#define __GOTTCHA_OPTIONS

#include <deque>
#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>

typedef enum {CIRCULAR, LINEAR, UNKNOWN_TOPOLOGY} Topology;

typedef enum {
	STRAIN_LEVEL,
	SPECIES_LEVEL,
	GENUS_LEVEL,
	FAMILY_LEVEL,
	ORDER_LEVEL,
	CLASS_LEVEL,
	PHYLUM_LEVEL,
	KINGDOM_LEVEL,
	NUM_TAXA_LEVEL
} TaxaLevel;

namespace std 
{
	// Allow the use of TaxaLevel as a hash key
	template <> struct hash<TaxaLevel>
	{
		size_t operator()(const TaxaLevel &m_taxa) const
		{
			return static_cast<std::size_t>(m_taxa);
		};
	};
}

struct Options
{
	public:
		unsigned int word_size; // k-mer signature size
		unsigned int min_output_frag; // Minimum genome fragement size
		bool compress_output;
		bool quit;
		Topology default_topology;
		std::string output_root_dir;
		std::unordered_map<TaxaLevel, std::string> taxa_level_output_prefix;
		
		enum {SILENT, INFORMATIVE, VERBOSE} verbose;

		std::unordered_map<TaxaLevel, std::string> taxa_level_mapping_file;
		std::deque<std::string> groups_to_squash;

		std::string log_file;

		// The maximum amount of RAM (in bytes) that each worker is assumed to possess
		size_t max_worker_RAM;
				
		Options(int argc, char* argv[], bool m_load_options);

		template<class T> friend 
			size_t mpi_size(const T &m_obj);
		template<class T> friend 
			unsigned char* mpi_pack(unsigned char* m_ptr, const T &m_obj);
		template<class T> friend 
			unsigned char* mpi_unpack(unsigned char* m_ptr, T &m_obj);
};

template<> size_t mpi_size(const Options &m_obj);
template<> unsigned char* mpi_pack(unsigned char* m_ptr, const Options &m_obj);
template<> unsigned char* mpi_unpack(unsigned char* m_ptr, Options &m_obj);

std::string taxa_level_name(const TaxaLevel &m_level);

#endif // __GOTTCHA_OPTIONS
