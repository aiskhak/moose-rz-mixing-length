#!/bin/bash
set -e

APP=/homes/aiskhak/projects/fv_app/fv_app-opt

for level in n16 n32 n64 n128; do
  echo "Running $level"
  cd $level
  rm -f *.csv *.e *.log
  $APP -i mms_rz_viscosity_${level}.i > run.log 2>&1
  tail -15 run.log
  cd ..
done

python3 check_convergence.py
