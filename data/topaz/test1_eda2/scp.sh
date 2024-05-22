#!/bin/bash

rsync -avP topaz:/group/director2183/data/test/eda2/20190913_051147/antenna_positions/* .
rsync -avP topaz: /group/director2183/data/test/eda2/20190913_051147/antenna_positions/20191104_033537_eda2_ch1_ant256_midday_avg1_vis_????_channel000_time000000.fits .
rsync -avP topaz:/group/director2183/data/test/eda2/20190913_051147/antenna_positions/doit! .

