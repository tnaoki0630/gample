defaults
==============

FlagForEquationDefault
----------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``Particle``
     - ``@0``
   * - ``EMField``
     - ``@0``
   * - ``MCCollision``
     - ``@0``

ParamForTimeIntegrationDefault
------------------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``Start``
     - ``@0``
   * - ``End``
     - ``@0``
   * - ``TimeStep``
     - ``@0``
   * - ``LogOutput``
     - ``@0``
   * - ``ParticleOutput``
     - ``@0``
   * - ``FieldOutput``
     - ``@0``
   * - ``ProgressOutput``
     - ``@0``
   * - ``ProjectName``
     - ``@"undefined"``
   * - ``RestartName``
     - ``@"undefined"``

ParamForComputingDefault
------------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``MaximumParticleNumber``
     - ``@0.08``
   * - ``ThreadGroupSize``
     - ``@256``
   * - ``IntegrationChunkSize``
     - ``@256``
   * - ``AggregationThreshold``
     - ``@0.1``
   * - ``ChebyshevDegree``
     - ``@5``
   * - ``AMGCycleType``
     - ``@1``
   * - ``MaxiterForPoisson``
     - ``@200``
   * - ``ToleranceForPoisson``
     - ``@(1e-7)``

ParamForParticleDefault
-----------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``ParticleName``
     - ``@"undefined"``
   * - ``pNumSetMethod``
     - ``@"undefined"``
   * - ``pNumValue``
     - ``@0``
   * - ``Charge``
     - ``@0``
   * - ``Mass``
     - ``@(-1)``
   * - ``WeightSetMethod``
     - ``@"undefined"``
   * - ``WeightValue``
     - ``@0``
   * - ``GenerateType``
     - ``@"undefined"``
   * - ``InitialPosX``
     - ``@[@(-1), @(-1)]``
   * - ``InitialPosY``
     - ``@[@(-1), @(-1)]``
   * - ``InitialVel``
     - ``@[@0, @0, @0]``
   * - ``InitialTemp``
     - ``@0``

BoundaryCoditionDefault
-----------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``RegionName``
     - ``@"undefined"``
   * - ``BCType``
     - ``@"undefined"``
   * - ``Value``
     - ``@0``

SourceForParticleDefault
------------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``ParticleName``
     - ``@"undefined"``
   * - ``GenerateType``
     - ``@"undefined"``
   * - ``SrcSetMethod``
     - ``@"undefined"``
   * - ``SrcVal``
     - ``@0``
   * - ``GeneratePosX``
     - ``@[@(-1), @(-1)]``
   * - ``GeneratePosY``
     - ``@[@(-1), @(-1)]``
   * - ``GenerateVel``
     - ``@[@0, @0, @0]``
   * - ``GenerateTemp``
     - ``@0``
   * - ``GenVelForElectron``
     - ``@[@0, @0, @0]``
   * - ``GenTempForElectron``
     - ``@0``

ParamForFieldDefault
--------------------

.. list-table::
   :header-rows: 1

   * - Key
     - Default
   * - ``NumberOfGridX``
     - ``@0``
   * - ``NumberOfGridY``
     - ``@0``
   * - ``GridSizeOfX``
     - ``@0``
   * - ``GridSizeOfY``
     - ``@0``
   * - ``InitializeTypeOfE``
     - ``@"undefined"``
   * - ``AmplitudeOfE``
     - ``@[@0, @0, @0]``
   * - ``FilePathOfEx``
     - ``@"undefined"``
   * - ``FilePathOfEy``
     - ``@"undefined"``
   * - ``FilePathOfEz``
     - ``@"undefined"``
   * - ``InitializeTypeOfB``
     - ``@"undefined"``
   * - ``AmplitudeOfB``
     - ``@[@0, @0, @0]``
   * - ``FilePathOfBx``
     - ``@"undefined"``
   * - ``FilePathOfBy``
     - ``@"undefined"``
   * - ``FilePathOfBz``
     - ``@"undefined"``
   * - ``WeightingOrder``
     - ``@5``
   * - ``BorisOrder``
     - ``@2``
   * - ``PtclPosOffset``
     - ``@0``
