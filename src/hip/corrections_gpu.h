#ifndef __BLINK_CORRECTIONS_GPU_H__
#define __BLINK_CORRECTIONS_GPU_H__

#include <astroio.hpp>
#include <memory_buffer.hpp>

void apply_geometric_corrections_gpu(Visibilities &xcorr, float *w_gpu, MemoryBuffer<double>& frequencies);

void apply_cable_lengths_corrections_gpu(Visibilities &xcorr, MemoryBuffer<double>& cable_lengths, MemoryBuffer<double>& frequencies);

#endif