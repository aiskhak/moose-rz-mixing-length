#!/bin/bash
set -e

APP=/homes/aiskhak/projects/fv_app/fv_app-opt

run_case () {
  level=$1
  np=$2

  echo "Running $level on $np MPI ranks"
  cd $level
  rm -f *.csv *.e *.log

  mpiexec -n $np $APP -i velocity_mms_${level}.i > run_np${np}.log 2>&1

  tail -25 run_np${np}.log
  cd ..
}

run_case n16  2
run_case n32  4
run_case n64  8
run_case n128 8

python3 check_convergence.py
