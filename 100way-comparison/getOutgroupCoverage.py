#!/usr/bin/env python
# Args: species tree, cactus log
import sys
import re
import networkx as nx
from sonLib.nxnewick import NXNewick

if len(sys.argv) != 3:
    print "Usage: %s tree log" % sys.argv[0]
    sys.exit(1)

tree = sys.argv[1]
log = open(sys.argv[2])

tree = NXNewick().parseString(tree)
nxTree = nx.Graph(tree.nxDg)
dm = nx.algorithms.shortest_paths.weighted.all_pairs_dijkstra_path_length(nxTree)
nameToId = {}
for node in nx.nodes(nxTree):
    nameToId[tree.getName(node)] = node

firstOutgroups = {}

print 'ingroup\tcoverage\tfirstOutgroup\tfirstOutgroupDist\tnumOutgroupsTotal'

ingroupData = {}

for line in log:
    line = line.strip()
    m = re.search(r'Cumulative coverage of (.*) outgroups on ingroup (.+)\.fa.*: (.*)', line)
    if m:
        num = m.group(1)
        ingroup = m.group(2)
        coverage = m.group(3)
        firstOutgroup = firstOutgroups[ingroup]
        ingroupData[ingroup] = (ingroup, coverage, firstOutgroup, dm[nameToId[firstOutgroup]][nameToId[ingroup]], num)
    m = re.search(r'Coverage on (.*)\.fa.* from outgroup #1, (.*)\.fa.*', line)
    if m:
        ingroup = m.group(1)
        firstOutgroups[ingroup] = m.group(2)

for ingroup, coverage, firstOutgroup, dist, num in ingroupData.values():
    print '%s\t%s\t%s\t%s\t%s' % (ingroup, coverage, firstOutgroup, dist, num)

