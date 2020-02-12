#ifndef __TAXA_ID
#define __TAXA_ID

#include <stdint.h>
#include <string>

// Extend the NCBI taxa id (integer) to handle strain id as well.
// TaxaId: <32 bit NCBI taxa id>.<32 bit strain id>

typedef uint64_t TaxaId;

TaxaId str_to_taxa_id(const std::string &m_tid, bool m_ignore_strain = false);
std::string taxa_id_to_str(const TaxaId &m_tid, bool m_show_strain = true);

#endif // __TAXA_ID
