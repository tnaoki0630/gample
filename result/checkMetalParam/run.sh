#!/bin/bash
for ((i=6;i<9;i++))
do

for ((j=0;j<i;j++))
do

export OMP_NUM_THREADS=4
echo "job tgs$((2**$j))ics$((2**($i))) startted."
./PIC -i tgs$((2**$j))ics$((2**($i))).json
wait %1

done
done
