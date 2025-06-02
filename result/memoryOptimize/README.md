- 積分用のバッファを Particle クラス側で確保していたが、同時に使うことはないので Moment に集約。
- メモリは20%ほど削減できたが、若干動作が不安定になってるかもしれない。理由はよくわからない。
![physMem](https://github.com/tnaoki0630/gample/blob/master/result/memoryOptimize/physMem.png)
![timeICDE](https://github.com/tnaoki0630/gample/blob/master/result/memoryOptimize/timeICDE.png)