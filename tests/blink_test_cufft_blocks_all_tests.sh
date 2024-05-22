#!/bin/bash

# 15 streams no blocks :
../tests/blink_test_cufft_blocks_all.sh NSTREAMS15_NOBLOCKS

# 1 streams and blocks ON :
../tests/blink_test_cufft_blocks_all.sh NSTREAMS1_BLOCKS "-s 1 -B"

# 15 streams and blocks ON :
../tests/blink_test_cufft_blocks_all.sh NSTREAMS15_BLOCKS "-s 15 -B"



