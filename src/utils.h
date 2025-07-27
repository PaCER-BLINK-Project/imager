#ifndef __BLINK_IMAGER_UTILS_H__
#define __BLINK_IMAGER_UTILS_H__

#include <memory_buffer.hpp>
#include <complex>
#include "images.hpp"

Images image_averaging_cpu(const Images& images);


template <typename T>
inline void compare_buffers(MemoryBuffer<T>& a, MemoryBuffer<T>&b ){
   for(size_t i {0}; i < a.size(); i++){
      if(std::abs(a[i] - b[i]) >= 1e-4) {
         std::cout << "Elements differ at position " << i << ": a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
         throw std::exception{};
      }
   }
}

template <>
inline void compare_buffers(MemoryBuffer<std::complex<float>>& a, MemoryBuffer<std::complex<float>>&b ){
   for(size_t i {0}; i < a.size(); i++){
      if(std::abs(a[i].real() - b[i].real()) >= 1e-4 || std::abs(a[i].imag() - b[i].imag()) >= 1e-4) {
         std::cout << "Elements differ at position " << i << ": a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
         //throw std::exception{};
      }
   }
}

void memdump(char *ptr, size_t nbytes, std::string filename);
#endif 