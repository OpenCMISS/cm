#!/usr/bin/env python

"""
Combine exnode files exported by OpenCMISS running in parallel into one exnode file.
This will read all "*.part*.exnode" files in the current directory and write out
a single "<name>.exnode" file where <name> is the name of the first part.exnode file found.

This script doesn't do anything about exelem files.

If you just want to visualise the solution in cmgui then read in all the the separate
exnode and exelem files. See examples/ClassicalField/Laplace/Laplace/visualse.com for example.
"""

import os

exnodes = [f for f in os.listdir(".") if f.find(".part") > -1 and f.find(".exnode") > -1]

if len(exnodes) == 0:
    raise RuntimeError, "No exnode files found"

print "Merging: "
for f in exnodes:
    print f,
print "\nOutput file: "
exnode_out_name = exnodes[0].split(".")[0]+".exnode"
print exnode_out_name
exnode_out = open(exnode_out_name,'w')

# get data from all exnode files
first=True
exnode_data=[]
for filename in exnodes:
    f = open(filename,'r')
    lines = f.readlines()
    if first:
        for (i,line) in enumerate(lines):
            if line.find("Node:") > -1:
                first_node = i
                break
        info = lines[0:first_node]
        exnode_out.writelines(info)
        first=False
    data = {}
    start = True
    for (i,line) in enumerate(lines):
        if line.find("Node:") > -1:
            node_num = int(line.split()[-1])
            data[node_num] = []
            if start:
                start=False
        if not start:
            data[node_num].append(line)
    exnode_data.append(data)
    f.close()

# build combined exnode file
all_nodes = {}
error_count = 0
for data in exnode_data:
    for k in data.keys():
        if k in all_nodes.keys():
            error_count += 1
            if error_count == 20:
                print "20 nodes already in data, not showing any more errors."
            elif error_count < 20:
                print "node", k, "is already in data"
    all_nodes.update(data)
nodes = all_nodes.keys()
nodes.sort()
for node in nodes:
    exnode_out.writelines(all_nodes[node])
