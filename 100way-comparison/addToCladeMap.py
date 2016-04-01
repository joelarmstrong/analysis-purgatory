#!/usr/bin/env python
# One-off to append to a clade-map file: a (genomeName, clade) matrix
# Args: halFile cladeName cladeMap
from sonLib.bioio import popenCatch, system
import sys

halFile = sys.argv[1]
cladeName = sys.argv[2]
cladeMapPath = sys.argv[3]

# first, read the clade map so we don't have duplicate names
alreadyMapped = set()
system("touch %s" % cladeMapPath)
cladeMap = open(cladeMapPath, "r")
for line in cladeMap:
    fields = line.split()
    alreadyMapped.add(fields[0])

cladeMap = open(cladeMapPath, "a")

names = popenCatch("halStats --genomes %s" % (halFile))
for name in names.split():
    if name not in alreadyMapped:
        cladeMap.write("%s\t%s\n" % (name, cladeName))
