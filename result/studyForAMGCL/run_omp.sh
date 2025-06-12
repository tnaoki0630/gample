#!/bin/bash
for ((i=0;i<6;i++))
do

echo "OMP=$(($i+1)) start"
export OMP_NUM_THREADS=$(($i+1))
./PIC -i data/condition.json
mv CG_log.xml CG_OMP$(($i+1))_log.xml
wait %1

done
