#include <astroio.hpp>
#include <memory_buffer.hpp>
#include "corrections.h"
#include "pacer_imager_defs.h"
#include <cmath>



void apply_geometric_corrections_cpu(Visibilities &xcorr, const MemoryBuffer<float>& w, const MemoryBuffer<double>& frequencies){
    if(xcorr.on_gpu()) xcorr.to_cpu();
    int n_ant = xcorr.obsInfo.nAntennas;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++){
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++) {
            const unsigned int n_baselines = static_cast<unsigned int>( (n_ant * (n_ant + 1)) / 2 );
            for(unsigned int baseline = 0; baseline < n_baselines; baseline++){
                unsigned int ant1 {static_cast<unsigned int>(-0.5 + std::sqrt(0.25 + 2*baseline))};
                unsigned int ant2 {baseline - ((ant1 + 1) * ant1)/2};
                double w_val = w[baseline];
                double angle = 2.0 * M_PI * w_val * frequencies[fine_channel] / SPEED_OF_LIGHT;
                double sin_angle, cos_angle;
                sincos(angle, &sin_angle, &cos_angle);
                for(int pol_idx {0}; pol_idx < 4; pol_idx++){ // TODO: avoid hardcoding 4 here..
                    std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2) + pol_idx;

                    double re = vis[0].real();
                    double im = vis[0].imag();

                    double re_prim = re * cos_angle - im * sin_angle;
                    double im_prim = im * cos_angle + re * sin_angle;

                    std::complex<double> vis_new(re_prim, im_prim);
                    *vis = vis_new;
                }
            }
        }
    }
}


void apply_cable_lengths_corrections_cpu(Visibilities &xcorr, const MemoryBuffer<double>& cable_lengths, const MemoryBuffer<double>& frequencies)
{
    if(xcorr.on_gpu()) xcorr.to_cpu();
    int n_ant = xcorr.obsInfo.nAntennas;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int time_step = 0; time_step < xcorr.integration_intervals(); time_step++)
    {
        for (int fine_channel = 0; fine_channel < xcorr.nFrequencies; fine_channel++)
        {
            for (int ant1 = 0; ant1 < n_ant; ant1++)
            {
                for (int ant2 = 0; ant2 <= ant1; ant2++)
                {
                    double cableDeltaLen = cable_lengths[ant2] - cable_lengths[ant1];
                    double angle = -2.0 * M_PI * cableDeltaLen * frequencies[fine_channel] / SPEED_OF_LIGHT;

                    double sin_angle, cos_angle;
                    sincos(angle, &sin_angle, &cos_angle);
                    for(int pol_idx {0}; pol_idx < 4; pol_idx++){ // TODO: avoid hardcoding 4 here..
                        std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2) + pol_idx;

                        double re = vis[0].real();
                        double im = vis[0].imag();
                        double re_prim = re * cos_angle - im * sin_angle;
                        double im_prim = im * cos_angle + re * sin_angle;

                        std::complex<double> vis_new(re_prim, im_prim);
                        *vis = vis_new;
                    }
                }
            }
        }
    }
}