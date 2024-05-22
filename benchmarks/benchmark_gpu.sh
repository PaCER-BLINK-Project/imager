#!/bin/bash

# PARAMETERS :
n_iterations=1000 # number of iterations (to make it run longer - maybe not the best way) - -P option
if [[ -n "$1" && "$1" != "-" ]]; then
   n_iterations=$1
fi

echo "##########################################"
echo "PARAMETERS:"
echo "##########################################"
echo "n_iterations = $n_iterations"
echo "##########################################"


# -P 1000 or so ...
module load blink_test_data/devel

cp ${BLINK_TEST_DATADIR}/eda2/20200209/chan_204_20200209T034646_vis_????.fits .
cp ${BLINK_TEST_DATADIR}/eda2/20200209/antenna_locations.txt .
cp ${BLINK_TEST_DATADIR}/eda2/20200209/calsol_merged.txt .

# 
pwd
which pacer_dirty_imager
# -P NUMBER_OF_LOOP ITERATIONS 
echo "time ./pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 1024 -w N -o miriad -c calsol_merged.txt -Z -v -100 -P ${n_iterations}"
time ./pacer_dirty_imager chan_204_20200209T034646 -f 159.375 -a antenna_locations.txt -n 1024 -w N -o miriad -c calsol_merged.txt -Z -v -100 -P ${n_iterations}
