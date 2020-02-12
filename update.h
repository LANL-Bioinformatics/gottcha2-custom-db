#ifndef __UPDATE
#define __UPDATE

#include <string>
#include <sstream>
#include <fstream>

class UpdateInfo : public std::stringstream {
	
	private:
		size_t buffer_size;
		std::ofstream flog;
		
	public:
	
		UpdateInfo(const std::string &m_log_file)
		{
			init( std::string() );
			
			Log(m_log_file);
		};
		
		UpdateInfo()
		{
			init( std::string() );
		};
		
		void Log(const std::string &m_log_file)
		{
			// Only attempt to open a log if the user has provided
			// a non-empty string
			if( !m_log_file.empty() ){
			
				flog.open( m_log_file.c_str() );

				if(!flog){
					throw __FILE__ ":Log: Unable to open the log file";
				}
			}
		};
		
		void init(const std::string &m_prefix);
		void flush();
		void close();
};

#endif // __UPDATE
