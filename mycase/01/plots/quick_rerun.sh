#!/bin/bash
set -e

echo "Rerunning the simulation."
sleep 5
echo "Removing old output files"
rm -rf plots/*
rm -rf *.out*
./cas_opt 01cas.in || true
echo "Visualize the results"
python3 plotting.py
echo "Finished"