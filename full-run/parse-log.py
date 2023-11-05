import sys
import pprint
from collections import defaultdict

filename = "full-run/full-run-{}.log".format(sys.argv[1])
count_dict = defaultdict(list)
count = 0
region = 1
sim = 1

with open(filename, "r") as f:
    for line in f:
        if "NegLL, Grad:" in line:
            last_value = line.strip().split(" ")[-1]
            if "e+" in last_value:
                count += 1
        elif "Finished region" in line:
            if count > 0:
                count_dict[count].append((sim, region))
            count = 0
            region += 1
        elif "Finished sim" in line:
            sim += 1
            region = 1


pprint.pprint(count_dict)