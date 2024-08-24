#ifndef __COMMON_H__
#define __COMMON_H__

#include <exception>
#include <string>
#include <sstream>

#define ENV_DATA_ROOT_DIR "BLINK_TEST_DATADIR"

class TestFailed : public std::exception {
    private:
    std::string message;

    public:
    TestFailed(const char *msg)  noexcept : std::exception(),  message {msg} {}
    virtual ~TestFailed() {};
    virtual const char* what() const noexcept override {return message.c_str();}
};


template <typename T>
bool complex_vectors_equal(const std::complex<T>* a, const std::complex<T>* b, size_t length){
    double delta;
    const double TOL {1e-5};
    for(size_t i {0}; i < length; i++){
        if (std::abs(a[i]) == 0) 
            delta = std::abs(b[i]);
        else 
            delta = std::abs(a[i] - b[i]);
        
        if (delta > TOL) {
            std::cout << "Elements at position " << i << " differs (delta = " << delta <<"): " << "a[i] = " << a[i] << ", b[i] = " << b[i] << std::endl;
            return false;
        }
    }
    return true;
}


void compare_xcorr_to_fits_file(Visibilities& xcorr, std::string filename){
    auto vis2 = Visibilities::from_fits_file(filename, xcorr.obsInfo);
    for(size_t int_time {0}; int_time < xcorr.integration_intervals(); int_time++){
        for(size_t fine_channel {0}; fine_channel < xcorr.nFrequencies; fine_channel++){
            size_t n_nans {0};
            size_t total {0};
            for(size_t a1 {0}; a1 < xcorr.obsInfo.nAntennas; a1++){
                for(size_t a2 {0}; a2 < a1; a2++){
                    std::complex<float> *p1 = xcorr.at(int_time, fine_channel, a1, a2);
                    std::complex<float> *p2 = vis2.at(int_time, fine_channel, a1, a2);
                    for(size_t p {0}; p < 4; p++){
                        total++;
                        if(std::isnan(p1->real()) && std::isnan(p2->real()) && std::isnan(p2->imag()) && std::isnan(p1->imag())){
                            n_nans++;
                            continue;
                        }
                        if(*p1 != *p2){
                            std::stringstream ss;
                            ss << "xcorr differs from " << filename << "!!!!" << std::endl;
                            ss << "[a1 = " << a1 << ", a2 = " << a2 << "] p1 = " << *p1 << ", p2 = " << *p2 << std::endl;
                            throw TestFailed(ss.str().c_str());
                        }
                    }
                }
            }
        }
    }
}


void load_dump(std::string filename, char *& buffer, size_t& size){
    std::ifstream infile (filename, std::ifstream::binary);
    // get size of file
    infile.seekg (0,infile.end);
    size = infile.tellg();
    infile.seekg (0);

    // allocate memory for file content
    buffer = new char[size];

    // read content of infile
    infile.read (buffer, size);
    infile.close();
}

#endif
