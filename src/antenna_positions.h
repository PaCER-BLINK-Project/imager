#ifndef _ANTENNA_POSITIONS__
#define _ANTENNA_POSITIONS__

#include <vector>
#include <string>
using namespace std;

// InputMapping in msfitslib now 
#include <array_config_common.h>
#include <memory_buffer.hpp>

class CObsMetadata;

/*class InputMapping
{
public :
   int input;
   int antenna;
   string szAntName;
   char pol;
   int delta;
   int flag;  
   
   double x;
   double y;
   double z;

   InputMapping() 
   : input(-1), antenna(-1), pol('U'), delta(0), flag(0), x(0), y(0), z(0)
   {};

//   static int read_mapping_file( std::vector<InputMapping>& inputs , const char* filename="instr_config.txt" );
};*/

class CAntennaPositions : public std::vector<InputMapping> 
{
public :
   CAntennaPositions( CObsMetadata* pMetaData=NULL );
   ~CAntennaPositions();
   
   int ReadAntennaPositions(const char* filename, bool bConvertToXYZ=false ); // default MRO location
   
   // if pMetaData=NULL    -> meta data not provided -> calculate UVW assuming zenith pointing (i.e. SKA-Low station)
   // else pMetaData!=NULL -> meta data provided -> calculate UVW using these metadata (e.g. for the MWA)
   int CalculateUVW(MemoryBuffer<float>& u_cpu, MemoryBuffer<float>& v_cpu, MemoryBuffer<float>& w_cpu, bool bIncludeAutos=false );
   
   // metadata information :
   CObsMetadata* m_pMetaData;
   
   // location of the observatory is part of the antenna location and relevant :
   static double m_ObservatoryLongitudeDeg;
   static double m_ObservatoryLattitudeDeg;
   static double m_ObservatoryHeight;
};

#endif 