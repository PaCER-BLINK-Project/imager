#include "pacer_imager.h"
#include "pacer_imager_thread.h"
// #include <mydate.h>

CImagingRequest::CImagingRequest()
: m_DataReal(NULL), m_DataImag(NULL), m_Ants(-1), m_Pols(-1)
{
}

CImagingRequest::CImagingRequest( float* data_real, 
                    float* data_imag,
                    int n_ant, 
                    int n_pol,
                    double frequency_mhz, 
                    int n_pixels,
                    double FOV_degrees,
                    double min_uv,      // minimum UV
                    bool do_gridding,    // excute gridding  (?)
                    bool do_dirty_image, // form dirty image (?)
                    const char* weighting, // weighting : U for uniform (others not implemented)
                    const char* szBaseOutFitsName
                  )
: m_DataReal(data_real), m_DataImag(data_imag), m_Ants(n_ant), m_Pols(n_pol), m_FrequencyMHz(frequency_mhz), m_Pixels(n_pixels), m_FOV_degrees(FOV_degrees), 
  m_MinUV(min_uv), m_DoGridding(do_gridding), m_DoDirtyImage(do_dirty_image), m_Weighting(weighting), m_szBaseOutFitsName(szBaseOutFitsName)
{

}


void CImagingQueue::AddImagingRequest( CImagingRequest& imaging_request )
{
   m_QueueLock.Lock();
   
   push( imaging_request );
   
   m_QueueLock.UnLock();   
   
   m_DataSignalPipe.Release();
}                 



CImagingRequest CImagingQueue::GetNextImagingReqeust()
{
   CImagingRequest ret;
   ret.m_DataReal = NULL;
   ret.m_DataImag = NULL;
   
   m_DataSignalPipe.Wait();
   
   m_QueueLock.Lock();
   
   if( size() <= 0 ){
      m_QueueLock.UnLock();
      return ret;
   }
      
   // get the next imaging request and remove from the front of the queue
   ret = front();
   pop();
   
   m_QueueLock.UnLock();
   
   return ret;
}

CImagingParameters::CImagingParameters( const char* szAntennaLocationFile, CImagingQueue* pImagingQueue )
: m_pImagingQueue(pImagingQueue)
{
   if( szAntennaLocationFile && strlen(szAntennaLocationFile) > 0 ){
      m_ImagerParameters.m_AntennaPositionsFile = szAntennaLocationFile;
   }      
}

CImagingParameters::CImagingParameters()
{
}

static void* RunImagingThread( void* pPtr )
{
   // CImagingQueue* pImagingQueue = (CImagingQueue*)pPtr;
   CImagingParameters* pImagingParameters = (CImagingParameters*)pPtr;
   if( !pImagingParameters ){
      printf("WARNING : pointer to imaging parameters has not been provided -> cannot execute the imaging thread\n");
      return NULL;
   }  
   CImagingQueue* pImagingQueue = pImagingParameters->m_pImagingQueue;
   
   if( !pImagingQueue ){
      printf("WARNING : pointer to imaging queue not provided -> cannot execute the imaging thread\n");
      return NULL;
   }
   
   bool bContinue=true;
   printf("INFO : starting imaging thread at unixtime = %d\n",(int)get_dttm());
   
   CPacerImager pacer_imager;
   pacer_imager.m_ImagerParameters = pImagingParameters->m_ImagerParameters;
//   if( strlen(pImagingParameters->m_ImagerParameters.m_AntennaPositionsFile.c_str()) > 0 ){
//      pacer_imager.m_ImagerParameters.SetParameters( pImagingParameters->m_ImagerParameters.m_AntennaPositionsFile.c_str() );
//   }
   pacer_imager.Initialise();
   
   while( bContinue ){
      // TODO :
      // add waiting for request to come (use pipe)
   
      CImagingRequest image_request = pImagingQueue->GetNextImagingReqeust();
      
      if( image_request.m_DataReal && image_request.m_DataImag ){
         printf("INFO : new imaging processing request received and being processed\n");
         
         // call imager here :
         bool imager_ret = pacer_imager.run_imager( image_request.m_DataReal, image_request.m_DataImag,
                                                    image_request.m_Ants, image_request.m_Pols,
                                                    image_request.m_FrequencyMHz, image_request.m_Pixels, 
                                                    image_request.m_FOV_degrees, image_request.m_MinUV,
                                                    image_request.m_DoGridding, image_request.m_DoDirtyImage,
                                                    image_request.m_Weighting.c_str(), image_request.m_szBaseOutFitsName.c_str() 
                                                  );
                                                  
         printf("DEBUG : pacer_imager.run_imager returned %d\n",imager_ret);
         
         // this is the only place where pointer to data should be removed after the imaging is done:
         delete [] image_request.m_DataReal;
         delete [] image_request.m_DataImag;
         
//      }else{
         // TODO : remove this and add proper waiting for request inside function GetNextImagingReqeust !!!
//         WaitMili(1);
      }
   }
   
   return (NULL);
}


CImagingThread::CImagingThread()
{}

CImagingThread::~CImagingThread()
{}

CImagingQueue& CImagingThread::GetImagingQueue()
{ 
   return m_ImagingQueue;
} 

bool CImagingThread::StartImagingThread( CImagingParameters* pImagingParameters )
{
   if( pImagingParameters ){
      if( !pImagingParameters->m_pImagingQueue ){
         pImagingParameters->m_pImagingQueue = &m_ImagingQueue;
      }
   }
   
   id = pthread_create( &m_ThreadID, NULL, RunImagingThread, (void*)pImagingParameters );

   printf("INFO : started imaging thread with id = %d, thread_id = %d\n",id,(int)m_ThreadID);   

   return true;
}   
