#!/bin/bash
set -e

echo "Testing automatic recompilation of the project when source files are modified."
echo "Run from the example directory, example (Aleiferis injector) used to set up the program."
echo "Removing old output files"
rm -rf plots/*
rm -rf *.out*
rm -rf cas_opt
echo "Removing old binary files"
rm -rf ../../bin/cas_opt
echo "New compilation: starting"
cd ../../src
make clean
make all_opt
echo "New compilation: done"
cd ../example/Al_inj 2>/dev/null || cd ../example/"Aleiferis injector"
cp ../../bin/cas_opt .
echo "Running the program again"
./cas_opt cas.in || true
echo "Visualize the results"
python3 result.py
echo "Rerun of the code is complete. Check the results."

