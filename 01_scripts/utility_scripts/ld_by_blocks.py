#!/usr/bin/env python3
"""Extract LD infos by blocks of size block_size

Usage:
    <program> input_ld block_size_in_kbp output_infos
"""

# Modules
from collections import defaultdict
import sys

# Parsing user input
try:
    input_ld = sys.argv[1]
    block_size = int(sys.argv[2]) * 1000
    output_infos = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Functions
def extract_summary(measures, value):
    infos = []
    measures_sorted = sorted(measures, reverse=True)

    positions = [0]
    one_hundred = 100
    spacing = (len(measures) - 1) / float(one_hundred)
    current_pos = 0
    for i in range(one_hundred):
        current_pos += spacing
        positions.append(int(current_pos))

    for i, pos in enumerate(positions):
        num = "000" + str(i * 10)
        num = num[-3:]
        infos.append((value + "_percent" + num, measures_sorted[pos], str(len(measures))))

    #for i in range(min(len(measures), 10)):
    #    num = "000" + str(i+1)
    #    num = num[-2:]
    #    infos.append((value + "_pos" + num, measures_sorted[i]))

    #infos.append(("num_measures", str(len(measures))))

    return infos

# Read LD file
first_line = True
pairs = defaultdict(lambda: defaultdict(list))
infos = defaultdict(lambda: defaultdict(list))

with open(input_ld) as infile:
    for line in infile:
        l = line.strip().split()

        if first_line:
            values = l[3:7]
            first_line = False

        else:
            p1 = int(l[0].split(":")[1]) // block_size 
            p2 = int(l[1].split(":")[1]) // block_size 
            positions = (p1, p2)

            for i, value in enumerate(l[3:7]):
                pairs[positions][values[i]].append(value)

# Extract summary infos
for p in pairs:
    for v in values:
        infos[p][v] = extract_summary(pairs[p][v], v)

# Write summary to output file
with open(output_infos, "wt") as outfile:
    for p in sorted(infos):
        for v in sorted(infos[p]):
            for m in infos[p][v]:
                stat, value, num_measures = m
                outfile.write("\t".join([str(p[0]), str(p[1]), stat, value, num_measures]) + "\n")
