#include <math.h>

#include "antenna_positions.h"
#include "pacer_geometry.h"
#include <myparser.h>
#include <myfile.h>
#include <mystrtable.h>
#include <bg_fits.h>
#include "observation_metadata.h"
#include <libnova_interface.h>

double CAntennaPositions::m_ObservatoryLongitudeDeg = MWA_LONGITUDE;
double CAntennaPositions::m_ObservatoryLattitudeDeg = MWA_LATTITUDE;
double CAntennaPositions::m_ObservatoryHeight = MWA_HEIGHT;


CAntennaPositions::CAntennaPositions( CObsMetadata* pMetaData )
: m_pMetaData(pMetaData)
{}

CAntennaPositions::~CAntennaPositions()
{}

int CAntennaPositions::CalculateUVW( CBgFits& fits_vis_u, CBgFits& fits_vis_v, CBgFits& fits_vis_w, 
                                     bool bSaveFits /*=false*/ , 
                                     const char* szOutDir /*="./"*/, 
                                     bool bIncludeAutos /* = false */ )
//                                     CObsMetadata* pMetaData /* = NULL */ )
{
   if( size() <= 0 ){
      return 0;
   }
   
   int n_ant = size();

   // initialise UVW arrays according to number of antennas :
   fits_vis_u.Realloc( n_ant , n_ant );
   fits_vis_v.Realloc( n_ant , n_ant );
   fits_vis_w.Realloc( n_ant , n_ant );

   // initalise :
   // set NaN to be consistent with CASA dump (was 0 but CASA dump script sets NaN for flagged antennas and upper half when option --conjugate is not used):
   fits_vis_u.SetNaN();
   fits_vis_v.SetNaN();
   fits_vis_w.SetNaN();


   double uxtime = 0.0; // TODO fix this hardcoding 
   double jd;
   double lmst_local = 0.00;
   PacerGeometry::UVWTimestepInfo uvwInfo; 
   if( m_pMetaData && m_pMetaData->m_bHasMetaFits ){
      uxtime = m_pMetaData->startUnixTime; // was 1419609944.00
      lmst_local = get_local_sidereal_time( uxtime, m_pMetaData->geo_long, jd )*15.00; // 116.670815 is geographic longigute in degrees
      // double get_local_sidereal_time(double uxtime_d,double geo_long_deg,double& jd_out)
      printf("DEBUG : CalculateUVW HA = %.6f [deg] , DEC = %.6f [rad] , lmst_local = %.8f [deg]\n",(m_pMetaData->haHrs*15.00),m_pMetaData->tilePointingDecRad,lmst_local);
      
      PacerGeometry::PrepareTimestepUVW( uvwInfo , *m_pMetaData );
   }
   int n_baselines = 0;
   
   for(int i=0;i<n_ant;i++){
      InputMapping& ant1 = (*this)[i];
      
      int start_j = (i+1);
      if( bIncludeAutos || true ){ // always calculate UVW for all including AUTOS 
         start_j = i;
      }
      
      for(int j=start_j;j<n_ant;j++){
         InputMapping& ant2 = (*this)[j];         
         
         double u = 0.00, v = 0.00, w = 0.00;
         if( m_pMetaData && m_pMetaData->m_bHasMetaFits ){
            double u2=0.00,v2=0.00,w2=0.00;
            double ha_rad = (m_pMetaData->haHrs*15.00)*(M_PI/180.00); // was - but testing with + was 313.37812500
//            double ha_rad = (313.37812500)*(M_PI/180.00);
            double dec_rad = m_pMetaData->tilePointingDecRad; // is -0.205727 radians 
            double lmst = lmst_local*(M_PI/180.00); // was 93.17821757
            double lmst2000 = lmst;

//            PacerGeometry::CalcUVW( ha_rad, lmst, dec_rad, lmst2000, ant1.x, ant1.y, ant1.z, u, v, w );
//            PacerGeometry::CalcUVW( ha_rad, lmst, dec_rad, lmst2000, ant2.x, ant2.y, ant2.z, u2, v2, w2 );
            // TEST:
            
            // TODO : uvwInfo.MJD + integrationtime/2 !!! 
            PacerGeometry::CalcUVW_NEW( uvwInfo , ant1.x, ant1.y, ant1.z, u, v, w );
            PacerGeometry::CalcUVW_NEW( uvwInfo , ant2.x, ant2.y, ant2.z, u2, v2, w2 );
           
            u = u - u2;
            v = v - v2;
            w = w - w2;
            PRINTF_DEBUG("DEBUG UVW : %d - %d : (%.3f,%.3f,%.3f) - (%.3f,%.3f,%.3f) = (%.4f,%.4f,%.4f)\n",i,j,ant1.x, ant1.y, ant1.z, ant2.x, ant2.y, ant2.z, u, v, w );
         }else{
            u = (ant1.x - ant2.x); // *(-1); TODO : image flip
            v = (ant1.y - ant2.y);
            w = (ant1.z - ant2.z);
         }
         
         // j,i (instead of i,j) to be consistent with CASA UVW array:
         fits_vis_u.setXY( j, i, u );
         fits_vis_v.setXY( j, i, v );
         fits_vis_w.setXY( j, i, w );
         
         // fill the other half of the array with -UVW values:
         fits_vis_u.setXY( i, j, -u );
         fits_vis_v.setXY( i, j, -v );
         fits_vis_w.setXY( i, j, -w );
         
         n_baselines++;
      }
   }
   
   if( bSaveFits ){
      char szFullPath[128];
      sprintf(szFullPath,"%s/u.fits",szOutDir);
   
      if( fits_vis_u.WriteFits( szFullPath ) ){
        printf("ERROR : could not write output file u.fits\n");        
      }
      
      sprintf(szFullPath,"%s/v.fits",szOutDir);
      if( fits_vis_v.WriteFits( szFullPath ) ){
        printf("ERROR : could not write output file v.fits\n");        
      }
      
      sprintf(szFullPath,"%s/w.fits",szOutDir);
      if( fits_vis_w.WriteFits( szFullPath ) ){
        printf("ERROR : could not write output file w.fits\n");        
      }
   }
   
   return n_baselines;
}

int CAntennaPositions::ReadAntennaPositions(const char* filename, bool bConvertToXYZ /*=false*/ )
{
   if( filename && strlen(filename) ){
      if( !MyFile::DoesFileExist( filename ) ){
         printf("ERROR: filename %s does not exist\n",filename);
         return -1;
      }
   }else{
      printf("ERROR : empty filename provided -> cannot continue\n");
      return -1;      
   }

   std::vector<InputMapping> tmp_antennas;

   clear();
   int n_ant_index = 0;
   MyFile file(filename);
   const char* pLine;
   if(!file.IsOpened()){
      file.Open( filename );
   }

   double geo_lat_rad = CAntennaPositions::m_ObservatoryLattitudeDeg*(M_PI/180.00);
   double geo_long_rad = CAntennaPositions::m_ObservatoryLongitudeDeg*(M_PI/180.00);

   bool remapping = false;
   while( (pLine = file.GetLine(TRUE)) ){
      if(strlen(pLine) && pLine[0]!='#'){
         MyParser pars=pLine;
         CMyStrTable items;
         pars.GetItems( items );
         
         if( strcmp(items[0].c_str(),"#") && items.size()>=4 ){
            InputMapping input;            
            input.szAntName = items[0].c_str();
            input.x = atof( items[1].c_str() );
            input.y = atof( items[2].c_str() );
            input.z = atof( items[3].c_str() );
            input.antenna = n_ant_index;
            if( items.size() >= 5 ){
               input.antenna = atof( items[4].c_str() );
               remapping = true;
            }
            
            input.local_x = input.x;
            input.local_y = input.y;
            input.local_z = input.z;
            
            if( CImagerParameters::m_bCalcEarthXYZ ){
               // for the MWA we require conversion from local (x,y,z) to Earth's (X,Y,Z) coordinates , for EDA2 we can keep using local (x,y,z)
               // TODO : in the future revise if this parameter is required or we can tell based on the data and other parameter (-Z option or -X ???)
               // see : /home/msok/bighorns/software/analysis/scripts/python/ant_local2earth.py
               double X0_casa = -2.55945e6; // hardcoded MWA in Earth XYZ coordinates, TODO : calculate using Geodetic2XYZ
               double Y0_casa = 5.09537e6;  // hardcoded MWA in Earth XYZ coordinates, TODO : calculate using Geodetic2XYZ
               double Z0_casa = -2.84906e6; // hardcoded MWA in Earth XYZ coordinates, TODO : calculate using Geodetic2XYZ
               
               double X0,Y0,Z0;
               PacerGeometry::Geodetic2XYZ( geo_lat_rad , geo_long_rad, CAntennaPositions::m_ObservatoryHeight, X0,Y0,Z0 );
               
               // see : /home/msok/bighorns/software/analysis/scripts/python/ant_local2earth.py
               input.x = X0 - input.local_y*sin(geo_lat_rad)*cos(geo_long_rad) - input.local_x*sin(geo_long_rad) + input.local_z*cos(geo_lat_rad)*cos(geo_long_rad);
               input.y = Y0 - input.local_y*sin(geo_lat_rad)*sin(geo_long_rad) + input.local_x*cos(geo_long_rad) + input.local_z*cos(geo_lat_rad)*sin(geo_long_rad);
               input.z = Z0 + input.local_y*cos(geo_lat_rad) + input.local_z*sin(geo_lat_rad);
               printf("DEBUG : %s Local (%.4f,%.4f,%.4f) -> Earth (%.4f,%.4f,%.4f) ( bConvertToXYZ = %d , X0,Y0,Z0 CASA (%.8f,%.8f,%.8f) vs. MINE (%.8f,%.8f,%.8f))\n",input.szAntName.c_str(),input.local_x,input.local_y,input.local_z,input.x,input.y,input.z,bConvertToXYZ,X0_casa,Y0_casa,Z0_casa,X0,Y0,Z0);
            }else{
               if( bConvertToXYZ ){ // add && false if CASA antenna coordinates are used as they are already in XYZ coordinate system 
                  PacerGeometry::ENH2XYZ_local( input.local_x, input.local_y, input.local_z, m_ObservatoryLattitudeDeg*(M_PI/180.00), input.x, input.y, input.z );                
                  printf("INFO : calculating XYZ in the coordinate system in the equatorial plane (%.4f,%.4f,%.4f) -> (%.4f,%.4f,%.4f)\n",input.local_x, input.local_y, input.local_z, input.x, input.y, input.z );
               }
            }
            
            tmp_antennas.push_back( input );
         }
      }
   }
   file.Close();        
   
   // 
   if( remapping ){
      // remapping :
      printf("INFO : antenna remapping is required (not 1-1)\n");
      resize( tmp_antennas.size() );
      for(int i=0;i<tmp_antennas.size();i++){
         (*this)[tmp_antennas[i].antenna] = tmp_antennas[i];
      }
   }else{
      // not remapping 
      printf("INFO : antenna list file is 1-1 -> no remapping is required\n");
      for(int i=0;i<tmp_antennas.size();i++){
         push_back( tmp_antennas[i] );
      }
   }
   
   if( true ){
      FILE* out_f = fopen("xyz.txt","w");
      for(int i=0;i<size();i++){
         fprintf(out_f,"%.8f %.8f %.8f %.8f %.8f %.8f\n",(*this)[i].x,(*this)[i].y,(*this)[i].z,(*this)[i].local_x,(*this)[i].local_y,(*this)[i].local_z);
      }
      fclose(out_f);
   }
   printf("INFO : read antenna positions from text file %s\n",filename);

   return size();
}


