#!/usr/bin/env python3
"""Extract LD infos by blocks of size block_size

Usage:
    <program> input_ld block_size_in_kbp output_infos
"""

# Modules
from collections import defaultdict
import gzip
import sys

# Defining functions
def myopen(_file, mode="rt"):
    if _file.endswith(".gz"):
        return gzip.open(_file, mode=mode)

    else:
        return open(_file, mode=mode)
# Parsing user input
try:
    input_ld = sys.argv[1]
    block_size = int(sys.argv[2]) * 1000
    output_infos = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read LD file
first_line = True
pairs = defaultdict(lambda: defaultdict(int))

# Gather R2 info by region blocks
with myopen(input_ld) as infile:
    for line in infile:
        l = line.strip().split()

        if first_line:
            values = l[6]
            first_line = False

        else:
            p1 = int(l[0].split(":")[1]) // block_size 
            p2 = int(l[1].split(":")[1]) // block_size
            positions = (p1, p2)
            rounded = round(float(l[6]), 3)
            pairs[positions][rounded] += 1

# Write percentile R2 values per block
with myopen(output_infos, "wt") as outfile:
    outfile.write("Pos1\tPos2\tPercent\tR2\n")

    for p in pairs:
        num_values = sum(pairs[p].values())
        p1 = p[0] * block_size // 1000 + block_size // 2000
        p2 = p[1] * block_size // 1000 + block_size // 2000

        for percentile in (x/200 for x in range(1, 11)):
            cumulative = 0.0

            for k in sorted(pairs[p], reverse=True):
                cumulative += pairs[p][k]

                if (cumulative / num_values) >= percentile:
                    outfile.write(f"{p1}\t{p2}\t{percentile:.4f}\t{k}\n")
                    break
