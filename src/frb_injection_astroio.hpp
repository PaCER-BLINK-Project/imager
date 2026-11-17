#pragma once
#include "frb_injection.hpp"
#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <vector>
namespace blink_frb {
// After phase/cable corrections, before gridding. Mutates science buffer in place.
// Return value: count of active time/channel slices, not number of detections.
std::size_t inject(Visibilities& vis, const Catalogue& catalogue,
                   double phase_ra_deg, double phase_dec_deg,
                   const MemoryBuffer<float>& u, const MemoryBuffer<float>& v,
                   const MemoryBuffer<float>& w, const MemoryBuffer<double>& frequencies,
                   const std::vector<int>& flags);
}
