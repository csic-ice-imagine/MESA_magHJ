#!/bin/bash

# Run with bash file.sh or ./file.sh (the latter with the right permissions).
# Change permissions with chmod +x file.sh to add execute permission.

echo "Cleaning up the mess..."

rm -f inlist_*MJ*
rm -f logfile
rm -f LOGS/*
rm -rf LOGS_*
rm -rf .mesa_temp_cache
rm -rf out
rm -f photos/*
rm -f planet_*
rm -rf __pycache__
rm -f testhub.yml
rm -f star

# echo "Done, bye bye!"
