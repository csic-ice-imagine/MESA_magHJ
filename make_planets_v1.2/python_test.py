#!/usr/bin/env python
import numpy as np

# Read test.
# Read as raw data.
a = np.genfromtxt("filename.txt")
# Specify names for columns (user-specified).
b = np.genfromtxt("filename.txt",names=["1","2","3","4","5","6","7","8","9"])
# Read names for columns (from the first commented line in the file).
c = np.genfromtxt("filename.txt",names=True)

# Dimensions.
print(a.shape)
print(b.shape)
print(c.shape)
print()

# Print some output.
print(a[0,0])
print(a[0,1])

print(b["1"][0])
print(b["2"][0])

print(c["AGE"][0])
print(c["RADIUS"][0])

print()
print("All done...")
