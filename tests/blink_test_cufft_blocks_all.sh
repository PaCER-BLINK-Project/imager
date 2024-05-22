#!/bin/bash

outdir="test"
if [[ -n "$1" && "$1" != "-" ]]; then
   outdir="$1"
fi

extra_options=""
if [[ -n "$2" && "$2" != "-" ]]; then
   extra_options=$2
fi

do_pretest=0
if [[ -n "$3" && "$3" != "-" ]]; then
   do_pretest=$3
fi


if [[ $do_pretest -gt 0 ]]; then
   # pre-test to check if images are correct (check re_00.fits and im_00.fits)
   echo "../tests/blink_test_cufft_blocks.sh 1 180 1 \"${extra_options}\""
   ../tests/blink_test_cufft_blocks.sh 1 180 1 "${extra_options}"
fi   

../benchmarks/benchmark_multi_buffer_full.sh 1 "180"  "1 2 5 10 20 30 32 40 50 100 200 300 400 500 1000 2000 3000 4000 5000" "${extra_options}"
../benchmarks/benchmark_multi_buffer_full.sh 1 "1024" "1 2 5 10 20 30 32 40 50 100 200 300 400 500 1000 2000 3000 4000 5000" "${extra_options}"
../benchmarks/benchmark_multi_buffer_full.sh 1 "4096" "1 2 5 10 20 30 32 40 50 100 200" "${extra_options}"



mkdir -p ${outdir}
echo "mv 0*txt ${outdir}/"
mv 0*txt ${outdir}/
