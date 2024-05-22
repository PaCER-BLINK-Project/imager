#!/bin/bash

n_iterations=1
if [[ -n "$1" && "$1" != "-" ]]; then
   n_iterations=$1
fi

image_size_list="180 1024 4096"
if [[ -n "$2" && "$2" != "-" ]]; then
   image_size_list="$2"
fi

n_images_list="1 2 5 10 20 30 40 50 100 200 300 400 500 1000 2000 3000 4000 5000"
if [[ -n "$3" && "$3" != "-" ]]; then
   n_images_list="$3"
fi

extra_options=""
if [[ -n "$4" && "$4" != "-" ]]; then
   extra_options="$4"
fi

for image_size in `echo ${image_size_list}`
do
   for n_images in `echo ${n_images_list}`
   do
      ok=1
      if [[ $image_size -ge 1024 ]]; then
         if [[ $n_images -gt 5000 ]]; then
            ok=0
         fi
      fi      
      
      if [[ $ok -gt 0 ]]; then
         outfile_name=`echo "$n_images $image_size $n_iterations" | awk '{printf("%05d_%05d_%05d.txt",$1,$2,$3);}'`
      
         echo "../benchmarks/benchmark_multi_buffer.sh $n_images $image_size $n_iterations \"${extra_options}\" > ${outfile_name} 2>&1"
         ../benchmarks/benchmark_multi_buffer.sh $n_images $image_size $n_iterations "${extra_options}" > ${outfile_name} 2>&1
      else
         echo "WARNING : n_images = $n_images and image_size = $image_size is too much -> skipped"
      fi
   done    
done
