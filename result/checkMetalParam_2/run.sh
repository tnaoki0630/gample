#!/bin/bash
for ((i=0;i<6;i++))
do

echo "tgs$((1*2**$i)) ics64 start"
./PIC -i result/checkMetalParam_2/condition_$((1*2**$i))_64.json -o result/checkMetalParam_2/log_$((1*2**$i))_64.xml &
wait %1

echo "tgs$((1*2**$i)) ics32 start"
./PIC -i result/checkMetalParam_2/condition_$((1*2**$i))_32.json -o result/checkMetalParam_2/log_$((1*2**$i))_32.xml &
wait %1

done
