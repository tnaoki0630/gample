## how to run
- run.sh を gample ディレクトリで実行
- wait %1 で同一ターミナルで投げられたジョブ番号（左の場合、1番上のバックグラウンドタスク）が終わるまで待機できる。
## memo
- ICDEの計算速度: tgs=1, ics=128 が良さそう <- reduction における待ち合わせがなくなる影響かも？ tgs = 4 とかで reduction をなくした場合にどうなるかみてもいいかもしれない。
    - threadGroupSize=1 が最速
    ![comp_ics64](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/comp_ics64.png)
    ![comp_ics128](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/comp_ics128.png)
    - integrateChunkSize=128 が最速
    ![comp_tgs1](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/comp_tgs1.png)
    ![comp_tgs2](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/comp_tgs2.png)
<!-- - update, solvePoisson の計算速度: あんま変わらなさそう
    - update
    ![update_ics128](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/update_ics128.png)
    ![update_tgs2](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/update_tgs2.png)
    - solvePoisson
    ![solvePoisson_ics128](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/solvePoisson_ics128.png)
    ![solvePoisson_tgs2](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/solvePoisson_tgs2.png) -->
- physicalMemoryUsage: arrSize*ics/tgs のメモリを確保するので、必ずしも tgs=1 がいいとは限らない。
    - tgsが小さいほど大きい。
    ![physicalFootprint_ics128](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/physicalFootprint_ics128.png)
    - icsが大きいほど大きい。
    ![physicalFootprint_tgs1](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/physicalFootprint_tgs1.png)
<!-- - meanPotential: けっこう違う
![meanPotential_ics128](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/meanPotential_ics128.png)
![meanPotential_tgs1](https://github.com/tnaoki0630/gample/blob/master/result/checkMetalParam/meanPotential_tgs1.png) -->