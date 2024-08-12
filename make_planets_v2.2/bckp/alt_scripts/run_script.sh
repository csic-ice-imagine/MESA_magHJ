#!/bin/bash

echo "Running sequential script..."

#
echo "Starting job 1..."
cp -f input/inlist_create_m1 inlist_create
./rn > output.txt
cp -r LOGS LOGS_m1
rm output.txt
echo "Job 1 done."

#
echo "Starting job 2..."
cp -f input/inlist_create_m2 inlist_create
./rn > output.txt
cp -r LOGS LOGS_m2
rm output.txt
echo "Job 2 done."

echo "All done!"
