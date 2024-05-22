#include "apply_calibration.h"
#include "pacer_imager.h"
#include "pacer_imager_defs.h"
#include <myparser.h>
#include <myfile.h>
#include <mystrtable.h>


CCalSol::CCalSol( int ant, std::complex<double> _xx, std::complex<double> _xy, std::complex<double> _yx, std::complex<double> _yy )
// : xx(_xx), xy(_xy), yx(_yx), yy(_yy) 
: m_ant(ant)
{
   m_cal[0] = _xx;
   m_cal[1] = _xy;
   m_cal[2] = _yx;
   m_cal[3] = _yy;
}

CCalSols::CCalSols()
: m_freq_channel(-1), m_frequency_mhz(-1)
{   
}

CCalSols::CCalSols(const CCalSols& right )
{
   (*this) = right;
}

CCalSols& CCalSols::operator=( const CCalSols& right )
{
   clear();
   
   m_filename = right.m_filename;
   for(int i=0;i<right.size();i++){
      push_back( right[i] );
   }   
   
   return (*this);
}


// read using filename from member variable m_filename    
int CCalSols::read_calsolutions_text(  double sign_phase /* = -1 */ , bool bInverseAmp /*=false */, bool bForce /*=false*/ )
{
   if( strlen(m_filename.c_str()) ){
      return read_calsolutions_text( m_filename.c_str() , bForce, sign_phase );
   }
   
   printf("WARNING : no file name with calibration solutions specified (data will not be calibrated)\n");
   
   return 0;
}

// read using filename from member variable m_filename and chosing a different reading function according to file extension :
int CCalSols::read_calsolutions( double sign_phase /* = -1 */ , bool bInverseAmp /*=false*/, bool bForce /*=false*/ )
{
   if( strlen(m_filename.c_str()) ){
      if( strstr(m_filename.c_str(),".txt") ){
         return read_calsolutions_text( sign_phase, bInverseAmp, bForce );
      }else{
         printf("ERROR : format of calibration solution is not supported (will not use calibration from file %s)\n",m_filename.c_str());
         return 0;
      }
   }
   
   printf("WARNING : no file name with calibration solutions specified (data will not be calibrated)\n");
   return 0;
}


int CCalSols::read_calsolutions_text( const char* filename , bool bForce,  double sign_phase /* = -1 */ , bool bInverseAmp /*=false*/ )
{
   // testing and debugging :
   // WARNING : phase sign changed here but later in apply calibration, gains are used a (GAIN) NOT 1/GAIN, code : std::complex<double> z_cal = (g1)*z*(g2_c); ( bool CPacerImager::ApplySolutions )
   // WARNING : this may need some revision and cleaning up, but at least this combination seems to give correct flux scale of the Sun - as in the MIRIAD images !
   // double sign_phase = -1; 
   printf("DEBUG : in CCalSols::read_calsolutions_text(%s , %.1f, %d)\n",filename,sign_phase,bInverseAmp);

   if( size() > 0 ){
      if( !bForce ){
         printf("INFO : calibration solutions already read earlier and it's not required to re-read them again\n");
         return size();
      }
   }

   if( filename && strlen(filename) ){
      if( !MyFile::DoesFileExist( filename ) ){
         printf("ERROR: filename %s does not exist\n",filename);
         return -1;
      }
   }else{
      printf("ERROR : empty filename provided -> cannot continue\n");
      return -1;
   }

   clear();
   int n_ant_index = 0;
   MyFile file(filename);
   const char* pLine;
   if(!file.IsOpened()){
      file.Open( filename );
   }
   while( (pLine = file.GetLine(TRUE)) ){
      if( strlen(pLine) ){ // was && pLine[0]!='#' 
         MyParser pars=pLine;
         CMyStrTable items;
         pars.GetItems( items );

         if( strcmp(items[0].c_str(),"#") && items.size()>=4 ){
            int ant = atol( items[0].c_str() );
         
            double amp_xx = atof( items[1].c_str() ); // or 1/amp_xx
            if( bInverseAmp ){
               amp_xx = 1.00 / amp_xx;
            }
            double phase_deg_xx = sign_phase*atof( items[2].c_str() );
            double phase_rad_xx = phase_deg_xx*(M_PI/180.00);
            double xx_re = amp_xx*cos(phase_rad_xx);
            double xx_im = amp_xx*sin(phase_rad_xx);
            std::complex<double> xx( xx_re , xx_im );
            

            double amp_yy = atof( items[3].c_str() ); // or 1/amp_yy
            if( bInverseAmp ){
               amp_yy = 1.00 / amp_yy;
            }
            double phase_deg_yy = sign_phase*atof( items[4].c_str() );
            double phase_rad_yy = phase_deg_yy*(M_PI/180.00);
            double yy_re = amp_yy*cos(phase_rad_yy);
            double yy_im = amp_yy*sin(phase_rad_yy);
            std::complex<double> yy( yy_re , yy_im );
            
            std::complex<double> xy( 0, 0 );
            std::complex<double> yx( 0, 0 );
            
            push_back( CCalSol(ant, xx, xy, yx, yy) );
         }else{
            // parse frequency and channel from header :
            // # frequency channel 204 = 159.3750 MHz
            // # ANT AMP_X PHASE_X AMP_Y PHASE_Y AMP_XY PHASE_XY AMP_YX PHASE_YX
            // printf("DEBUG : skipping comment |%s|\n",pLine);
            if( items.size()>=7 ){
               int len = strlen("# frequency channel");
               if( strlen(pLine)>len && strncmp( pLine, "# frequency channel", len ) == 0 && items.size()>=7 ){
                  m_freq_channel = atol( items[3].c_str() );
                  m_frequency_mhz = atof( items[5].c_str() );
                  printf("DEBUG : frequency channel / value = %d / %.4f MHz\n",m_freq_channel,m_frequency_mhz);               
               }
            }
         }
      }
   }
   file.Close();
   printf("DEBUG : read %d solutions\n",int(size()));

   return size();
}

void CCalSols::show()
{
   if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL){   
      printf("\nCalibration solutions :\n");
      printf("ANT : CAL_xx_mag CAL_xx_phase_deg | CAL_yy_mag CAL_yy_phase_deg\n");
   
      for(int i=0;i<size();i++){
         CCalSol& calsol = (*this)[i];
      
         std::complex<double> cal_xx = calsol.m_cal[0];
         std::complex<double> cal_yy = calsol.m_cal[3];
      
         printf("%d : %.8f %.2f | %.8f %.2f\n",calsol.m_ant,std::abs(cal_xx),std::arg(cal_xx)*(180.00/M_PI),std::abs(cal_yy),std::arg(cal_yy)*(180.00/M_PI));
      }
   }
}
