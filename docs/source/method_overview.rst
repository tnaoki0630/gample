PIC-MCC overview
================

kinetic plasma
--------------
- ホールスラスタにおけるプラズマ異常輸送は運動論的な不安定性に由来している可能性がある．:cite:p:`Boeuf2017`

Particle-in-cell method
-----------------------
.. mermaid::

    flowchart TD
        init["初期値 $$ \mathbf x_{s,k}, \mathbf v_{s,k} $$"] --> loopStart[*]
        loopStart --> updateFld["電場更新 $$\nabla^2 \phi(\mathbf x) = -\frac{\rho(\mathbf x)}{\varepsilon_0} $$"]
        updateFld --> outputFld["output $$ \rho, \mathbf E, \mathbf B $$"]
        updateFld --> updatePtcl["粒子位置更新 $$ \mathbf v_{s,k} += \frac{q_{s,k}}{m_{s,k}}(\mathbf E(\mathbf x_{s,k}) + \mathbf v_{s,k}\times\mathbf B(\mathbf x_{s,k})) $$"]
        updatePtcl --> outputMoment["output $$ n_s, \mathbf u_s, \mathbf P_s $$"]
        updatePtcl --> deposit["電荷密度更新 $$\rho(\mathbf x) = \sum_s\sum_k q_{s,k}w_{s,k}\delta(\mathbf x - \mathbf x_p) $$"]

Monte_Carlo_Collision
---------------------
- piyo