#!/usr/bin/env python
# gets a TSV of sequence loss relative to its children (as ratio of
# avg genome size) for each ancestor.
import sys
from sonLib.bioio import popenCatch
from sonLib.nxnewick import NXNewick
import networkx as nx

if len(sys.argv) < 2:
    print "Usage: %s halFile" % sys.argv[0]
    sys.exit(1)

halFile = sys.argv[1]

stats = popenCatch("halStats %s" % halFile)

treeString = stats.split("\n")[2]
tree = NXNewick().parseString(treeString)
nxTree = nx.Graph(tree.nxDg)
dm = nx.algorithms.shortest_paths.weighted.all_pairs_dijkstra_path_length(nxTree)
nameToId = {}
for node in nx.nodes(nxTree):
    nameToId[tree.getName(node)] = node

genomeData = [row.split(", ") for row in stats.split("\n")[5:] if row != '']

lengths = {}

for row in genomeData:
    name = row[0]
    length = int(row[2])
    lengths[name] = length

print 'name\tsequenceLossRatio\tdistBtwnChildren'

for row in genomeData:
    name = row[0]
    numChildren = int(row[1])
    if numChildren == 0:
        continue
    children = [child.strip() for child in popenCatch("halStats --children %s %s" % (name, halFile)).split(" ")]
    childLengths = [lengths[child] for child in children]
    avgChildLength = float(sum(childLengths))/len(childLengths)
    myLength = lengths[name]
    sequenceLossRatio = 1 - myLength/avgChildLength
    assert len(children) == 2
    dist = dm[nameToId[children[0]]][nameToId[children[1]]]
    print "%s\t%s\t%s" % (name, sequenceLossRatio, dist)

