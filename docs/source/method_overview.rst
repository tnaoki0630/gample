PIC-MCC overview
================

kinetic plasma
--------------
- ホールスラスタにおけるプラズマ異常輸送は運動論的な不安定性に由来している可能性がある．:cite:p:`Boeuf2017`

Particle-in-cell method
-----------------------
.. mermaid::

   flowchart TD
     A["入力 $$x_0$$"] --> B{"更新 $$x_{k+1}=g(x_k)$$"}
     B -->|$$\lVert r_k\rVert <\varepsilon$$| C["収束 $$x^*$$"]
     B -->|else| B

Monte_Carlo_Collision
---------------------
- piyo