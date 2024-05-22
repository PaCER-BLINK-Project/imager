# Building the project in this branch is as follows:

  - Laptop/workstation using NVIDIA :
    - CPU :
      cd build
      cmake ..
      make 
      make test
      sudo make install

    - GPU :
      cd build_gpu
      cmake -DUSE_HIP=ON -DCUDA_SM=86 .. # make use to change CUDA_SM to value relevant to your specific GPU (check on https://en.wikipedia.org/wiki/CUDA )
      make
      make test
      sudo make install

  - Setonix - at least temporarily there are separate build_setonix.sh and CMakeLists.txt_SETONIX scripts for Setonix/HIP architecture
    - CPU version
      ./build.sh cpu - "-DENABLE_PROFILER=ON"

    - GPU version
      salloc -p gpu --gres=gpu:1 --time 0:30:00 -n1 --account=???
      ./build.sh gpu - "-DENABLE_PROFILER=ON  -DUSE_HIP=ON"


# pacer_dirty_imager

   Pacer BLINK imaging software (including test and demonstration programs). Fast imaging program developed in C++ for the BLINK PaCER project. 

# required 

   sudo apt-get install libnova-dev fftw2 fftw-dev fftw3-dev libhdf5-dev libcfitsio-dev

   msfitslib follow the installation instructions at https://github.com/marcinsokolowski/msfitslib 

# installation:

   On Topaz super computer :

   git clone git@146.118.67.64:blink/imager.git

   ./build.sh

On any desktop computer:

   mkdir build

   cd build  

   cmake ..
   ( for debug version : cmake -DCMAKE_BUILD_TYPE=Debug .. )

   make

   sudo make install

# current programs :

   - The package includes script to dump correlation matrix from CASA or UVFITS files

      - Starting with demonstration software using CASA dumps of visibility and UVW data
      - Scripts for dumping CASA ms are in cotter_wsclean/scripts/casa/


   - pacer_uvgrid_generator - visibility simulator (filling UV grid with requested pattern) - it's just for testing.


   - pacer_dirty_image - reads visibilities and UVW data from FITS files as produced by script cotter_wsclean/scripts/casa/get_corrmatrix_from_casa.py (may be moved here now)
     then gridds the data on UV grid, performs invFFT and gets dirty image. 
     Has already been demonstrated to work well with EDA2 data (SKA-Low station data) and currently work is underway to make it working with MWA data.

     Example usage :

           2049x2049 image : pacer_dirty_image_test 20191104_033537_eda2_ch1_ant256_midday_avg1 -p _channel000_time000000 -n 2049

     generate visibilties and do inverse FFT :

                pacer_uvgrid_generator -A 100.00 -n 129
                pacer_dirty_image in out -g 0 -r uv_grid_re.fits -i uv_grid_im.fits

     or 8193x8193 :
                pacer_uvgrid_generator -A 100.00 -n 8193
                pacer_dirty_image in out -g 0 -r uv_grid_re.fits -i uv_grid_im.fits -n 8193
                          

Test data:

   EDA2 test dataset for Topaz (taken from as in scp.sh script):

     cd data/topaz/test1_eda2

     ./doit_topaz!

     The program is expected to create :

        Image of the sky named : dirty_image_20220406T035747506.fits

        Imaginary image named : 20220406T035747506_imag.fits

     WARNING : the name conventions are being changed. So, they are likely to change soon, especially the name of imaginary part is just WRONG now !

     The reference images are included in the test dataset, hence 

        dirty_image_20220406T035747506.fits should be equal to dirty_image_20220406T035747506_real_TEMPLATE_512x512.fits

        20220406T035747506_imag.fits should be equal to dirty_image_20220406T035747506_imag_TEMPLATE_512x512.fits

   For any desktop computer just execute ./doit! script which should give the same results as above.

       
