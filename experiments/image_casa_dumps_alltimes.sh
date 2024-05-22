#!/bin/bash

for timestep in `ls -d ??????`
do
   cd ${timestep}
   image_casa_dumps.sh
   cd ..
done
