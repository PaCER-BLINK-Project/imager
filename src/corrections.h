#ifndef __BLINK_CORRECTIONS_H__
#define __BLINK_CORRECTIONS_H__

#include <astroio.hpp>
#include <memory_buffer.hpp>
void apply_geometric_corrections_cpu(Visibilities &xcorr, const MemoryBuffer<float>& fits_vis_w, const MemoryBuffer<double>& frequencies);


void apply_cable_lengths_corrections_cpu(Visibilities &xcorr, const MemoryBuffer<double>& cable_lengths, const MemoryBuffer<double>& frequencies);


#endif