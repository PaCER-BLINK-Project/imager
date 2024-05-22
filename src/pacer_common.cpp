#include "pacer_common.h"
#include <sys/time.h>
#include <stdio.h>

double get_dttm_decimal()
{
   time_t ut_time = 0;
   int usec = 0; 
   timeval t_info;
   
   gettimeofday( &t_info , NULL );
   // double ret = t_info.tv_sec + t_info.tv_usec/1000000.00;
   ut_time = (time_t)t_info.tv_sec;
   usec = t_info.tv_usec;
   
   double ret = ut_time + usec/1000000.00;
   return ret;
}

// additinal local functions:
const char* get_ut_string( double ux_time, char* out_buffer, const char* formater /*="%.2u%.2u%.2u_%.2u%.2u%.2u.%.3u"*/)
{
   time_t ut_time = 0;
   int usec = 0; 
   
   if( ux_time<=0.00001 ){ // just in case <= 0.00 does not work well with float
      timeval t_info;
      gettimeofday( &t_info , NULL );
      // double ret = t_info.tv_sec + t_info.tv_usec/1000000.00;
      ut_time = (time_t)t_info.tv_sec;
      usec = t_info.tv_usec;
   }else{
      ut_time = int(ux_time);
      usec = int( (ux_time - double(ut_time))*1000000.00 );
   }
   printf("DEBUG : get_ut_string uxtime = %.6f -> ut_time = %d, usec = %d\n",ux_time,int(ut_time),usec);

   struct tm gmtime_tm;
   if(gmtime_r( &ut_time, &gmtime_tm )){
      // bug ??? first %.2u -> %.4u ???      
      int millisec = int(usec/1000.00);
      sprintf(out_buffer, formater, gmtime_tm.tm_year+1900,(gmtime_tm.tm_mon+1),gmtime_tm.tm_mday,gmtime_tm.tm_hour,gmtime_tm.tm_min,gmtime_tm.tm_sec,millisec);
   }
   
   return out_buffer;
}

const char* get_filename( double ux_time, char* out_buffer,
                          const char* full_dir_path /*="./"*/, const char* prefix /*="dirty_image"*/, 
                          const char* postfix /*=""*/, const char* formater /*="%.2u%.2u%.2uT%.2u%.2u%.2u.%.3u" */ )
{
   char date_time[64];

   get_ut_string( ux_time , date_time , formater );   
   sprintf(out_buffer,"%s/%s%s%s.fits",full_dir_path,prefix,date_time,postfix);
   
   return out_buffer;
}

const char* get_filename_base( double ux_time, char* out_buffer,
                          const char* full_dir_path /*="./"*/, const char* prefix /*="dirty_image"*/, 
                          const char* postfix /*=""*/, const char* formater /*="%.2u%.2u%.2uT%.2u%.2u%.2u.%.3u" */ )
{
   char date_time[64];

   get_ut_string( ux_time , date_time , formater );   
   sprintf(out_buffer,"%s/%s%s%s",full_dir_path,prefix,date_time,postfix);
   
   return out_buffer;
}

