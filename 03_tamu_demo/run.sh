#!/bin/bash
source ~/miniforge/bin/activate moose
stdbuf -oL -eL mpirun -np 4 /homes/aiskhak/projects/fv_app/fv_app-opt -i tamu.i > log.csv 2>&1