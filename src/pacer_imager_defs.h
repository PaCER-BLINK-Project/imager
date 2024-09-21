#ifndef _PACER_IMAGER_DEFS_H__
#define _PACER_IMAGER_DEFS_H__

// high-time resolution clock
#include <chrono>
using namespace std::chrono;

// visibility type :
#define VISIBILITY_TYPE float

// Physical constants :
#define VEL_LIGHT  299792458.0    //!< speed of light in vacuo. meters per second
#define SPEED_OF_LIGHT 299792458.0        // speed of light in m/s

#define MWA_LATTITUDE -26.703319        // Array latitude. degrees North
#define MWA_LONGITUDE 116.67081         // Array longitude. degrees East
#define MWA_HEIGHT 377               // Array altitude. meters above sea level


// Program constants :
# define MAX_VIS 1e12

// verbosity level :
#define IMAGER_DEBUG_LEVEL        5
#define IMAGER_INFO_LEVEL         3
#define IMAGER_BENCHMARK_LEVEL    2
#define IMAGER_WARNING_LEVEL      1
#define IMAGER_ERROR_LEVEL        0
#define IMAGER_ALL_MSG_LEVEL      100
#define IMAGER_DEEP_DEBUG_LEVEL   200 

// level of saving intermediate and test files
#define SAVE_FILES_NONE  0
#define SAVE_FILES_FINAL 1
#define SAVE_FILES_FINAL_ALL 2 // all including Imaginary final image
#define SAVE_FILES_INFO  3 
#define SAVE_FILES_DEBUG 4
#define SAVE_FILES_ALL   100

#define PRINTF_DEBUG   if(CPacerImager::m_ImagerDebugLevel>=IMAGER_DEBUG_LEVEL)printf
#define PRINTF_BENCHMARK    if(CPacerImager::m_ImagerDebugLevel>=IMAGER_BENCHMARK_LEVEL)printf
#define PRINTF_INFO    if(CPacerImager::m_ImagerDebugLevel>=IMAGER_INFO_LEVEL)printf
#define PRINTF_WARNING if(CPacerImager::m_ImagerDebugLevel>=IMAGER_WARNING_LEVEL)printf
#define PRINTF_ERROR   if(CPacerImager::m_ImagerDebugLevel>=IMAGER_ERROR_LEVEL)printf

#ifdef _PACER_PROFILER_ON_

// OLD using clock -> NEW using chrono

#ifdef _OBSOLATE_USE_CLOCK_
   #define PACER_PROFILER_START clock_t t1=clock(),t2;
   #define PACER_PROFILER_RESTART t1=clock();
   #define PACER_PROFILER_STOP clock_t t2=clock();
   #define PACER_PROFILER_STOP2 t2=clock();
   #define PACER_PROFILER_END(desc) t2=clock(); \
                                    mystring msg=get_clock_in_sec_string( t2-t1 );\
                                    printf("PROFILER %s %s\n",desc,msg.c_str());
   #define PACER_PROFILER_SHOW(desc) t2=clock(); \
                                     msg=get_clock_in_sec_string( t2-t1 );\
                                     printf("PROFILER %s %s\n",desc,msg.c_str());
#else                                  
   #define PACER_PROFILER_START_TIMER1   high_resolution_clock::time_point timer1_start = high_resolution_clock::now();\
                                         high_resolution_clock::time_point timer1_end;
                                         
   #define PACER_PROFILER_START   high_resolution_clock::time_point t1 = high_resolution_clock::now();
   #define PACER_PROFILER_RESTART t1=high_resolution_clock::now();
   #define PACER_PROFILER_STOP    high_resolution_clock::time_point t2 = high_resolution_clock::now();
   #define PACER_PROFILER_STOP2   t2 = high_resolution_clock::now();
   #define PACER_PROFILER_END(desc) high_resolution_clock::time_point t2 = high_resolution_clock::now();\
                                    duration<double> time_span_gridding = duration_cast<duration<double>>(t2 - t1);\
                                    mystring msg=get_chrono_in_sec_string( time_span_gridding.count() );\
                                    printf("PROFILER %s %s\n",desc,msg.c_str());
                                 
   #define PACER_PROFILER_SHOW(desc) t2 = high_resolution_clock::now();\
                                     time_span_gridding = duration_cast<duration<double>>(t2 - t1);\
                                     msg=get_chrono_in_sec_string( time_span_gridding.count() );\
                                     printf("PROFILER %s %s\n",desc,msg.c_str());                                                                   

   #define PACER_PROFILER_SHOW_TIMER1(desc) timer1_end = high_resolution_clock::now();\
                                     time_span_gridding = duration_cast<duration<double>>(timer1_end - timer1_start);\
                                     msg=get_chrono_in_sec_string( time_span_gridding.count() );\
                                     printf("PROFILER %s %s\n",desc,msg.c_str());                                                                   
#endif
                                  
#else
#define PACER_PROFILER_START_TIMER1
#define PACER_PROFILER_START
#define PACER_PROFILER_RESTART
#define PACER_PROFILER_END(desc)
#define PACER_PROFILER_STOP
#define PACER_PROFILER_STOP2
#define PACER_PROFILER_SHOW(desc)
#define PACER_PROFILER_SHOW_TIMER1(desc)
#endif
                                 
enum eTelescopName_Type { eUnknown=0, eMWA=1, eEDA2=2, eAAVS2=3 };

enum eConvolutionKernel { eNoConvKernel=0, eConvKernelGauss=1, eConvKernelSinc=2, eConvKernelSincGauss=3 };


#endif 
