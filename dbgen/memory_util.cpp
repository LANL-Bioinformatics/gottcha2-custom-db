#include "gottcha_db.h"

using namespace std;

// Global variables for MPI
extern int mpi_numtasks;
extern int mpi_rank;

// Please note that the available RAM calculation currently assumes that each compute node is running
// on a separate physical host. This is not necessarily the case. This function needs to be updated so
// that the RAM available to each host is computed and then shared amoung the compute ranks running 
// on the host...
size_t available_RAM()
{
	ifstream fin("/proc/meminfo");
	
	if(!fin){
		throw __FILE__ ":available_RAM: Unable to open /proc/meminfo for reading! Is this a linux-based system?";
	}
	
	string line;
	
	size_t ret = 0;
	
	while( getline(fin, line) ){
		
		if(line.find("MemAvailable:") == string::npos){
			continue;
		}
		
		stringstream ssin(line);
		
		// MemAvailable:
		if( !(ssin >> line) ){
			throw __FILE__ ":available_RAM: Unable to parse header";
		}
		
		size_t mem = 0;
		
		if( !(ssin >> mem) ){
			throw __FILE__ ":available_RAM: Unable to parse memory";
		}
		
		string units;
		
		if( !(ssin >> units) ){
			throw __FILE__ ":available_RAM: Unable to parse memory units";
		}
		
		if(units != "kB"){
			throw __FILE__ ":available_RAM: Units are not in kB -- the /proc/meminfo format has changed!";
		}
		
		// The current format for /proc/meminfo specifies the available RAM in units of kB. Convert this
		// into bytes
		ret = mem*1024;
		
		return ret;
	}
	
	throw __FILE__ ":available_RAM: Unable to parse /proc/meminfo\nOld kernels should specify max ram via --RAM";
	return 0;
}

// Check all ranks to see if they have room to store more data in RAM
bool room_for_more(const bool &m_has_room)
{
	int has_room = m_has_room;
	int all_have_room = 0;
	
	// Use MPI_MIN to identify when one or more ranks no longer has room!
	MPI_Allreduce(&has_room, &all_have_room, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
	
	return all_have_room;
}
