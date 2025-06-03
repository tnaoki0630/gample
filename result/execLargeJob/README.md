- 処理が軽くなったことをいいことに、initPPC を 50 にして、weight 軽めの重い計算を流した。
- 1sec/cycle ほどで順調に動いていたが、50000サイクルほど経過すると 100sec オーダーかかるようになった。
![elapsedTime](https://github.com/tnaoki0630/gample/blob/master/resut/execLargeJob/elapsedTime.png)
- メモリを調べると rss が 25000 サイクルあたりから非線形に下がっているところがあり、とても怪しい挙動。physMem と rss は傾向が一致するはず。
![elapsedTime](https://github.com/tnaoki0630/gample/blob/master/resut/execLargeJob/elapsedTime.png)
- ChatGPT に心当たりがないか聞いたところ、物理メモリが不足した際にメモリの圧縮がかかるらしい。
- ガクッと下がってるのはメモリが圧縮されたタイミングで、50000 サイクルあたりは完全に諦めてストレージを仮想メモリとして使ってる感じらしい。
- だいたい合計1100万粒子くらいが限界かも？