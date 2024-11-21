#include "images.hpp"
#include "files.hpp"

namespace {
        void save_fits_file(const std::string filename, float* data, long side_x, long side_y){
        FITS fitsImage;
        FITS::HDU hdu;
        hdu.set_image(data,  side_x, side_y);
        // hdu.add_keyword("TIME", static_cast<long>(obsInfo.startTime), "Unix time (seconds)");
        // hdu.add_keyword("MILLITIM", msElapsed, "Milliseconds since TIME");
        // hdu.add_keyword("INTTIME", integrationTime, "Integration time (s)");
        // hdu.add_keyword("COARSE_CHAN", obsInfo.coarseChannel, "Receiver Coarse Channel Number (only used in offline mode)");
        fitsImage.add_HDU(hdu);
        fitsImage.to_file(filename);
    }
}

void Images::to_fits_files(const std::string& directory_path) {
    if(on_gpu()) to_cpu();
    for(size_t interval {0}; interval < this->integration_intervals(); interval++){
        for(size_t fine_channel {0}; fine_channel < this->nFrequencies; fine_channel++){
            float *current_data {this->data() + this->image_size() * this->nFrequencies * interval + fine_channel * this->image_size()}; 
            std::stringstream full_directory;
            full_directory << directory_path << "/" << "start_time_" << obsInfo.startTime << \
                "/" << "int_" << interval << "/coarse_" << obsInfo.coarseChannel << "/fine_ch" << fine_channel;
            std::string full_directory_str = full_directory.str();
            blink::imager::create_directory(full_directory_str);
            ::save_fits_file(full_directory_str + "/image_real.fits", current_data, this->side_size, this->side_size);
        }
    }
}