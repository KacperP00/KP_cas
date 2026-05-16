#!/bin/bash
set -e

echo "Automatic recompilation of the project when source files are modified."
sleep 5
echo "Run from the example directory, example (Aleiferis injector) used to set up the program."
echo "Removing old output files"
rm -rf plots/*
rm -rf *.out*
rm -rf cas_opt
echo "Removing old binary files"
rm -rf ../../bin/cas_opt
echo "New compilation: starting"
sleep 5
cd ../../src
make clean
make all_opt
echo "New compilation: done"
cd ../mycase/01
cp ../../bin/cas_opt .
echo "Running the program again"
./cas_opt cas.in || true
echo "Visualize the results"
python3 plotting.py
echo "Rerun of the code is complete. Check the results."

