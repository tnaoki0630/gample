#!/bin/bash
for ((j=0;j<9;j++))
do

export OMP_NUM_THREADS=4
echo "job tgs$((2**$j)) startted."
./PIC -i tgs$((2**$j)).json
wait %1

done
