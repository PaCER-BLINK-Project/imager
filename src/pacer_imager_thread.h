#ifndef _PACER_IMAGER_THREAD_H__
#define _PACER_IMAGER_THREAD_H__

#include <queue>
#include <string>
#include <mylock.h>
#include <mypipe.h>

#include "pacer_imager.h"

using namespace std;

class CImagingRequest
{
public :
   CImagingRequest();

   CImagingRequest( float* data_real, 
                    float* data_imag,
                    int n_ant, 
                    int n_pol,
                    double frequency_mhz, 
                    int n_pixels,
                    double FOV_degrees,
                    double min_uv=-1000,      // minimum UV
                    bool do_gridding=true,    // excute gridding  (?)
                    bool do_dirty_image=true, // form dirty image (?)
                    const char* weighting="", // weighting : U for uniform (others not implemented)
                    const char* szBaseOutFitsName=NULL
                  );


   float* m_DataReal;
   float* m_DataImag;
   int    m_Ants;   
   int    m_Pols;
   double m_FrequencyMHz;
   int    m_Pixels;   
   double m_FOV_degrees;
   double m_MinUV;
   bool   m_DoGridding;
   bool   m_DoDirtyImage;
   string m_Weighting;
   string m_szBaseOutFitsName;
   
};

class CImagingQueue : public queue<CImagingRequest>
{
protected :
   CMyMutex  m_QueueLock;
   CPipeLock m_DataSignalPipe; // giving signal that data are there or not 

public :
   void AddImagingRequest( CImagingRequest& imaging_request );   
   CImagingRequest GetNextImagingReqeust();
};

class CImagingParameters
{
public : 
  CImagingParameters();
  CImagingParameters( const char* szAntennaLocationFile, CImagingQueue* pImagingQueue );
  
//  string m_szAntennaLocationFile;
  CImagerParameters m_ImagerParameters;
  CImagingQueue* m_pImagingQueue;
};

class CImagingThread 
{
protected :
   CImagingQueue m_ImagingQueue;
//   CPacerImager  m_PacerImager;

public :
   pthread_t m_ThreadID;
   int       id;
   
   CImagingThread();
   ~CImagingThread(); 

   CImagingQueue& GetImagingQueue();
   virtual bool StartImagingThread( CImagingParameters* pImagingParameters );
};


#endif
