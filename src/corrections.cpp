#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <bg_fits.h>
#include "corrections.h"
#include "pacer_imager_defs.h"


void apply_geometric_corrections_cpu(Visibilities &xcorr, CBgFits &fits_vis_w, const MemoryBuffer<double>& frequencies){
    if(xcorr.on_gpu()){
        std::cout << "xcorr is on GPU!!! Moving to cpu.." << std::endl;
        xcorr.to_cpu();
    }else{
        std::cout << "xcorr is on CPU.." << std::endl;
    }
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/01_before_geo_corrections.fits");
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
                    double w = fits_vis_w.getXY(ant1, ant2);
                    double angle = - 2.0 * M_PI * w * frequencies[fine_channel] / SPEED_OF_LIGHT;
                    double sin_angle, cos_angle;
                    sincos(angle, &sin_angle, &cos_angle);

                    std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2);

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
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/01_after_geo_corrections.fits");
}

void apply_cable_lengths_corrections_cpu(Visibilities &xcorr, const MemoryBuffer<double>& cable_lengths, const MemoryBuffer<double>& frequencies)
{
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/02_before_cable_corrections.fits");
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

                    std::complex<VISIBILITY_TYPE> *vis = xcorr.at(time_step, fine_channel, ant1, ant2);

                    double re = vis[0].real();
                    double im = vis[0].imag(); // TODO : why do I need 1 here ???
                    double re_prim = re * cos_angle - im * sin_angle;
                    double im_prim = im * cos_angle + re * sin_angle;

                    std::complex<double> vis_new(re_prim, im_prim);
                    *vis = vis_new;
                }
            }
        }
    }
    //::compare_xcorr_to_fits_file(xcorr, "/scratch/director2183/cdipietrantonio/1276619416_1276619418_images_cpu_reference_data/1592584200/133/000/02_after_cable_corrections.fits");
}