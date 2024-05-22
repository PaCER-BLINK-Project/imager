#!/bin/bash

#!/bin/bash

# DEBUG : /home/msok/Desktop/PAWSEY/PaCER/software/cotter/build_20230818_debug
# r 1276619416_20200619163000_gpubox01_00.fits 1276619416_20200619163000_gpubox02_00.fits 1276619416_20200619163000_gpubox03_00.fits 1276619416_20200619163000_gpubox04_00.fits 1276619416_20200619163000_gpubox05_00.fits 1276619416_20200619163000_gpubox06_00.fits 1276619416_20200619163000_gpubox07_00.fits 1276619416_20200619163000_gpubox08_00.fits 1276619416_20200619163000_gpubox09_00.fits 1276619416_20200619163000_gpubox10_00.fits 1276619416_20200619163000_gpubox11_00.fits 1276619416_20200619163000_gpubox12_00.fits 1276619416_20200619163000_gpubox13_00.fits 1276619416_20200619163000_gpubox14_00.fits 1276619416_20200619163000_gpubox15_00.fits 1276619416_20200619163000_gpubox16_00.fits 1276619416_20200619163000_gpubox17_00.fits 1276619416_20200619163000_gpubox18_00.fits 1276619416_20200619163000_gpubox19_00.fits 1276619416_20200619163000_gpubox20_00.fits 1276619416_20200619163000_gpubox21_00.fits 1276619416_20200619163000_gpubox22_00.fits 1276619416_20200619163000_gpubox23_00.fits 1276619416_20200619163000_gpubox24_00.fits 


# cotter -absmem 64 -j 12 -timeres 1 -freqres 0.01 -edgewidth 80 -noflagautos -flagantenna 25,58,71,80,81,92,101,108,114,119,125 -m 20200619163000.metafits -noflagmissings -allowmissing -offline-gpubox-format -initflag 0  -centre 18h33m41.89s -03d39m04.25s -o 1276619416_20200619163000.ms  -nocablelength -nogeom 1276619416_20200619163000*gpubox*.fits
# -usepcentre
#
# Without : -nogeom -nocablelength
/home/msok/Desktop/PAWSEY/PaCER/software/cotter/build_20230908/cotter -noalign -absmem 64 -j 12 -timeres 1 -freqres 0.01 -edgewidth 0 -noflagautos -m 20200619163000.metafits -noflagmissings -allowmissing -offline-gpubox-format -initflag 0 -o 1276619416_20200619163000.ms -norfi -noantennapruning -nosbgains -noflagdcchannels -sbpassband subband-passband-32ch-unitary.txt 1276619416_20200619163000*gpubox*.fits > cotter.out 2>&1

echo "casapy --nologger -c /home/msok/github/pacer/software/imager/scripts/get_corrmatrix_from_casa.py 1276619416_20200619163000.ms --data_column=DATA --save_full_matrix --channel=0"
casapy --nologger -c /home/msok/github/pacer/software/imager/scripts/get_corrmatrix_from_casa.py 1276619416_20200619163000.ms --data_column=DATA --save_full_matrix --channel=0
