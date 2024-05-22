#ifndef _PACER_COMMON_H__
#define _PACER_COMMON_H__

#include <time.h>

// simple functions required for some internal calculations : date/time/filename etc :

// for default files name :
const char* get_ut_string( double ux_time , char* out_buffer, const char* formater="%.2u%.2u%.2u_%.2u%.2u%.2u%.3u" );
const char* get_filename(  double ux_time , char* out_buffer, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix="", const char* formater="%.2u%.2u%.2uT%.2u%.2u%.2u%.3u" );
const char* get_filename_base(  double ux_time , char* out_buffer, const char* full_dir_path="./", const char* prefix="dirty_image_", const char* postfix="", const char* formater="%.2u%.2u%.2uT%.2u%.2u%.2u%.3u" );

// get unixtime with decimal value:
double get_dttm_decimal();

#endif
