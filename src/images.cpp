#include "images.hpp"
#include "files.hpp"
#include <iomanip>

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

void Images::to_fits_files(const std::string& directory_path, bool save_as_complex, bool save_imaginary) {
    if(on_gpu()) to_cpu();
    MemoryBuffer<float> img_real(this->image_size(), false, false);
    MemoryBuffer<float> img_imag(this->image_size(), false, false);
    blink::imager::create_directory(directory_path);
    for(size_t interval {0}; interval < this->integration_intervals(); interval++){
        for(size_t fine_channel {0}; fine_channel < this->nFrequencies; fine_channel++){
            std::complex<float> *current_data {this->data() + this->image_size() * this->nFrequencies * interval + fine_channel * this->image_size()}; 
            std::stringstream full_directory;
            full_directory << std::setfill('0') << directory_path << "/" << "start_time_" << obsInfo.startTime << \
                "_" << "int_" << std::setw(2) << interval <<  std::setw(0) << "_coarse_" <<  std::setw(3) <<  obsInfo.coarseChannel <<  std::setw(0) << "_fine_ch" <<  std::setw(2) << fine_channel;
            std::string full_directory_str {full_directory.str()};
            if(save_as_complex){
                std::string filename {full_directory_str + "_image.fits"};
                std::complex<float>* p_data = this->at(interval, fine_channel);
                ::save_fits_file(filename, reinterpret_cast<float*>(p_data), this->side_size, this->side_size *2);
            }else{
                for(size_t i {0}; i < this->image_size(); i++){
                    img_real[i] = current_data[i].real();
                }
                ::save_fits_file(full_directory_str + "_image_real.fits", img_real.data(), this->side_size, this->side_size);
                if(save_imaginary){
                    for(size_t i {0}; i < this->image_size(); i++){
                        img_imag[i] = current_data[i].imag();
                    }
                    ::save_fits_file(full_directory_str + "_image_imag.fits", img_imag.data(), this->side_size, this->side_size);
                }
            }
        }
    }
}