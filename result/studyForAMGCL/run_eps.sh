#!/bin/bash
for ((i=0;i<6;i++))
do

export OMP_NUM_THREADS=4

# echo "threshold = 0.01*$((2**$i)), cycle = 1, estimation = True"
# ./PIC_estTrue -i eps$((2**$i))cyc1.json
# mv CG_log.xml CG_estTrue_eps$((2**$i))cyc1_log.xml
# wait %1

# echo "threshold = 0.01*$((2**$i)), cycle = 2, estimation = True"
# ./PIC_estTrue -i eps$((2**$i))cyc2.json
# mv CG_log.xml CG_estTrue_eps$((2**$i))cyc2_log.xml
# wait %1

# echo "threshold = 0.01*$((2**$i)), cycle = 1, estimation = False"
# ./PIC_estFalse -i eps$((2**$i))cyc1.json
# mv CG_log.xml CG_estFalse_eps$((2**$i))cyc1_log.xml
# wait %1

# echo "threshold = 0.01*$((2**$i)), cycle = 2, estimation = False"
# ./PIC_estFalse -i eps$((2**$i))cyc2.json
# mv CG_log.xml CG_estFalse_eps$((2**$i))cyc2_log.xml
# wait %1

echo "threshold = 0.01*$((2**$i)), cycle = 1, estimation = True"
./PIC_aggrEmin -i eps$((2**$i))cyc1.json
mv CG_log.xml CG_aggrEmin_eps$((2**$i))cyc1_log.xml
wait %1

echo "threshold = 0.01*$((2**$i)), cycle = 2, estimation = True"
./PIC_aggrEmin -i eps$((2**$i))cyc2.json
mv CG_log.xml CG_aggrEmin_eps$((2**$i))cyc2_log.xml
wait %1

done
