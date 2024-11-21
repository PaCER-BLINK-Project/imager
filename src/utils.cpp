#include "images.hpp"


Images image_averaging_cpu(const Images& images){
    size_t image_size = images.image_size();
    size_t n_images = images.integration_intervals() * images.nFrequencies;
    MemoryBuffer<float> avg_image {image_size};
    memset(avg_image.data(), 0, sizeof(float) * image_size);
    for(size_t img_id {0}; img_id < n_images; img_id++){
        const float *img_data = images.data() + img_id * image_size; 
        for(size_t i {0}; i < image_size; i++) avg_image[i] += img_data[i];
    }
    for(size_t i {0}; i < image_size; i++) avg_image[i] /= static_cast<float>(n_images);
   return {std::move(avg_image), images.obsInfo, images.obsInfo.nTimesteps, images.obsInfo.nFrequencies, images.side_size};
}