#include "utils.h"


Images image_averaging_cpu(const Images& images){
    size_t image_size = images.image_size();
    size_t n_images = images.integration_intervals() * images.nFrequencies;
    MemoryBuffer<std::complex<float>> avg_image {image_size};
    memset(avg_image.data(), 0, sizeof(std::complex<float>) * image_size);
    for(size_t img_id {0}; img_id < n_images; img_id++){
        const std::complex<float> *img_data = images.data() + img_id * image_size; 
        for(size_t i {0}; i < image_size; i++){
            float val = img_data[i].real();
            if(i == 0) printf("img_id = %ld, val = %f\n", img_id, val);
            avg_image[i].real(avg_image[i].real() + val);
        }
    }
    for(size_t i {0}; i < image_size; i++){
        avg_image[i].real(avg_image[i].real() / static_cast<float>(n_images));
    }
    printf("avg = %f\n", avg_image[0].real());
   return {std::move(avg_image), images.obsInfo, images.obsInfo.nTimesteps, images.obsInfo.nFrequencies, images.side_size};
}


void memdump(char *ptr, size_t nbytes, std::string filename){
    std::ofstream outfile;
    outfile.open(filename, std::ofstream::binary);
    outfile.write(ptr, nbytes);
    if(!outfile){
        std::cerr << "Error while dumping data to " << filename << std::endl;
        exit(1);
    }
    outfile.close();
}
