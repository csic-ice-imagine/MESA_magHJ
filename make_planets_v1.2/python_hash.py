#!/usr/bin/env python
from pathlib import Path
import math
import os
import shutil

# Add hash signs to the beginnings of text lines (for plotting with gnuplot).

# Define the current path.
folder_path = os.getcwd()
print()
print("Checking in path:",folder_path)
print()
print("Files identified for copying:")

# Recursively and sequentially (in a sorted way) list all files in all folders.
# Locate the files that start with a given string.
for folder_path, currentDirectory, files in sorted(os.walk(folder_path)):
   for file in files:
      # if file.startswith("history.dat"):
      # if file == "history.data":
      if ".data" in file and "hash" not in file:
         # out_file = file.replace(".data","_hash.data")
         out_file = "hash_" + file
         print(file,out_file)

         # Full file name including path.
         in_file_path = os.path.join(folder_path,file)
         out_file_path = os.path.join(folder_path,out_file)

         # Open the input and output files.
         input_file = open(in_file_path,'r')
         output_file = open(out_file_path,'w')

         ll = 0
         for line in input_file:
            output_line = line
            # Edit the first 6 lines (adding #).
            ll = ll + 1
            # The fourth line is empty (only contains a newline string, \n).
            if ll == 4:
               output_line = '#' + output_line
               # print("boing",ll)
            # All other lines.
            elif ll < 7:
               output_line = '#' +  output_line[1:]
               # output_line = '#' +  output_line
               # print("dough",ll)
            output_file.write(output_line)

         input_file.close()
         output_file.close()

print()
print("All done...")
