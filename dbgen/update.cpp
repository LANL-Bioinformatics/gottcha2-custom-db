#include "update.h"
#include <iostream>

using namespace std;

void UpdateInfo::init(const string &m_prefix)
{
	cerr << m_prefix;
	
	buffer_size = 0;
}

void UpdateInfo::flush()
{
	// Clear the buffer
	for(size_t i = 0;i < buffer_size;i++){
		cerr << '\b';
	}
	
	for(size_t i = 0;i < buffer_size;i++){
		cerr << ' ';
	}
	
	for(size_t i = 0;i < buffer_size;i++){
		cerr << '\b';
	}
	
	const string tmp = str();
	
	buffer_size = tmp.size();
	
	cerr << tmp;
	
	if( flog.is_open() ){
		flog << tmp << endl;
	}
	
	// Clear the stringstream buffer
	str( string() );
}

void UpdateInfo::close()
{
	cerr << endl;
	
	// Clear the stringstream buffer
	str( string() );
}
