#pragma once
#include "frb_sim.hpp"
#include <string>
#include <vector>
#include <map>
#include <complex>
#include <cstddef>

namespace blink_frb {
// A catalogue entry uses topocentric UTC Unix seconds, not GPS or barycentric time.
struct Event {
    std::string observation_id, id;
    double ra_deg, dec_deg;
    long double arrival_unix_s;
    double dm, width_s, peak_jy, ref_hz, units_per_jy;
    int phase_sign; // -1: inverse-FFT convention; +1: Andy's exponent
    int tsub, fsub; // tsub=0: exact integration in time; otherwise midpoint sampling
    // Optional scattering parameters; peak_jy is the intrinsic boxcar peak.
    double scattering_tau_s = 0.0;
    double scattering_ref_hz = 150e6; // separate from dispersion ref_hz
    double scattering_index = -4.0;
};
using Catalogue = std::vector<Event>;
Catalogue load_catalogue(const std::string& path);
void validate(const Event& event);

struct Geometry {
    std::string observation_id;
    long double start_unix_s;
    double dt_s, df_hz, phase_ra_deg, phase_dec_deg;
    std::size_t nt;
    std::vector<double> frequencies_hz;
};
struct Contribution {
    frbsim::LMN lmn;
    double amplitude; // already in visibility units
    int phase_sign;
};
// Key = time * nChannels + channel; only active cells are stored.
using Plan = std::map<std::size_t, std::vector<Contribution>>;
Plan make_plan(const Catalogue&, const Geometry&);
// Add unpolarized FRB(s) to [baseline][XX,XY,YX,YY].
// Baselines flagged by the imager are untouched. Non-finite existing
// parallel hands remain untouched. No cross-hand or auto changes.
void add_slice(std::complex<float>* slice, std::size_t nb,
               const float* u, const float* v, const float* w,
               const int* flags, double frequency_hz,
               const std::vector<Contribution>& contributions);
} // namespace blink_frb
