#!/bin/bash
for ((i=1;i<9;i++))
do

for ((j=1;j<=i;j++))
do

export OMP_NUM_THREADS=4
echo "job tgs$((2**$j))ics$((2**($i+1))) startted."
./PIC -i tgs$((2**$j))ics$((2**($i+1))).json
wait %1

done
done
