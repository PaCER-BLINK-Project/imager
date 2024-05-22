#ifndef _PACER_IMAGER_HIP_DEFINES_H
#define _PACER_IMAGER_HIP_DEFINES_H

#include "../pacer_imager_defs.h"

#define NTHREADS 1024

// implemeneted in pacer_imager_hip.cu
// void __cuda_check_error(gpuError_t err, const char *file, int line);

//#define CUDA_CHECK_ERROR(X)({\
//        __cuda_check_error((X), __FILE__, __LINE__);\
//})

#endif