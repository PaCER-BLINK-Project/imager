#ifndef __BLINK_IMAGES_H__
#define __BLINK_IMAGES_H__

#include <astroio.hpp>
#include <memory_buffer.hpp>
#include <complex>
#include <string>

class Images : public MemoryBuffer<float> {
    public:
    ObservationInfo obsInfo;
    unsigned int nIntegrationSteps;
    unsigned int nAveragedChannels;
    unsigned int nFrequencies;
    unsigned int side_size;

   Images(MemoryBuffer<float>&& data, const ObservationInfo& obsInfo, unsigned int nIntegrationSteps,
            unsigned int nAveragedChannels, unsigned int side_size) : MemoryBuffer {std::move(data)} {
        this->obsInfo = obsInfo;
        this->nIntegrationSteps = nIntegrationSteps;
        this->nAveragedChannels = nAveragedChannels;
        this->nFrequencies = obsInfo.nFrequencies / nAveragedChannels;
        this->side_size = side_size;
    }

    Images(const Images& other) : MemoryBuffer {other} {
        obsInfo = other.obsInfo;
        nIntegrationSteps = other.nIntegrationSteps;
        nFrequencies = other.nFrequencies;
        nAveragedChannels = other.nAveragedChannels;
        side_size = other.side_size;
    }

    Images(Images&& other) : MemoryBuffer {std::move(other)} {
        obsInfo = other.obsInfo;
        nIntegrationSteps = other.nIntegrationSteps;
        nFrequencies = other.nFrequencies;
        nAveragedChannels = other.nAveragedChannels;
        side_size = other.side_size;
    }

    Images& operator=(Images& other){
        if(this == &other) return *this;
        MemoryBuffer::operator=(other);
        obsInfo = other.obsInfo;
        nIntegrationSteps = other.nIntegrationSteps;
        nFrequencies = other.nFrequencies;
        nAveragedChannels = other.nAveragedChannels;
        side_size = other.side_size;
        return *this;
    }

    Images& operator=(Images&& other){
        if(this == &other) return *this;
        obsInfo = other.obsInfo;
        nIntegrationSteps = other.nIntegrationSteps;
        nFrequencies = other.nFrequencies;
        nAveragedChannels = other.nAveragedChannels;
        side_size = other.side_size;
        MemoryBuffer::operator=(std::move(other));
        return *this;
    }


    float *at(unsigned int interval, unsigned int frequency) {
        const size_t nValuesInTimeInterval {image_size() * nFrequencies};
        float *pData = this->data() + nValuesInTimeInterval * interval + image_size() * frequency;
        return pData;
    }

    /**
     * Number of time intervals integrated over by the correlator.
     */
    size_t integration_intervals() const {
        return (obsInfo.nTimesteps + nIntegrationSteps - 1) / nIntegrationSteps;
    }
    
    // Number of pixels in a single image.
    size_t image_size() const {
       return side_size * side_size;
    }

    // Number of images `data` array.
    size_t size() const {
        return this->integration_intervals() * nFrequencies;
    }
    
    /**
     * @brief Save visibilities to a FITS file on disk.
     * 
     * @param filename name of the output file.
     */
    void to_fits_file(const std::string& filename);


   void to_fits_files(const std::string& directory_path);
};


#endif