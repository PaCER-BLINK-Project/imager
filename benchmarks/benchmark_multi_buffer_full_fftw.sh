#!/bin/bash

n_iterations=1
if [[ -n "$1" && "$1" != "-" ]]; then
   n_iterations=$1
fi

n_streams=1
if [[ -n "$2" && "$2" != "-" ]]; then
   n_streams=$2
fi

image_sizes="180 1024 4096"
if [[ -n "$3" && "$3" != "-" ]]; then
   image_sizes="$3"
fi


echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "n_iterations = $n_iterations"
echo "n_streams = $n_streams"
echo "image_sizes = $image_sizes"
echo "##########################################"

# 8192 -> 4096 to make it the same as GPU benchmarks 
for image_size in `echo ${image_sizes}`
do
   for n_images in `echo "1 10 20 30 40 50 100 200 300 400 500 1000 2000 3000 4000 5000"`
   do
      ok=1
      if [[ $image_size -ge 1024 ]]; then
         if [[ $n_images -gt 950 ]]; then
            ok=0
         fi
      fi      
      
      if [[ $ok -gt 0 ]]; then
         outfile_name=`echo "$n_images $image_size $n_iterations" | awk '{printf("%05d_%05d_%05d.txt",$1,$2,$3);}'`
      
         echo "../benchmarks/benchmark_multi_buffer_fftw.sh $n_images $image_size $n_iterations $n_streams > ${outfile_name} 2>&1"
         ../benchmarks/benchmark_multi_buffer_fftw.sh $n_images $image_size $n_iterations $n_streams > ${outfile_name} 2>&1
      else
         echo "WARNING : n_images = $n_images and image_size = $image_size is too much -> skipped"
      fi
   done    
done

