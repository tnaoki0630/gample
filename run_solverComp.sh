#!/bin/bash
solver="ICDbfrSP"
./PIC_$solver -i data/debye_$solver.json -o test
wait %1

solver="ICDaftSP"
./PIC_$solver -i data/debye_$solver.json -o test

