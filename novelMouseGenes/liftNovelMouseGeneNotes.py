#!/usr/bin/env python
import csv
import re
import sys
from collections import defaultdict
from glob import glob
clusters = []

with open('cluster.csv') as clustercsv:
    clusterReader = csv.reader(clustercsv)
    for row in clusterReader:
        cluster = []
        for entry in row:
            if entry == '':
                break
            fields = entry[1:-1].split(',')
            genome = fields[0]
            gene = fields[1]
            cluster.append((genome, gene))
        clusters.append(cluster)

# map from genome name -> gene names (ordered)
genes = defaultdict(list)
# map from genome name -> gene name -> annotation (which may be empty)
annotations = defaultdict(dict)

paths = glob('Downloads/mouse strains - *.novel+paralog.csv')
for path in paths:
    genome = re.match(r'Downloads/mouse strains - (.*).novel\+paralog.csv', path).group(1)
    with open(path) as f:
        reader = csv.reader(f)
        gene = None
        for i, row in enumerate(reader):
            empty = row == [] or row[0] == ''
            if i % 4 == 0:
                gene = row[0].split()[0]
                genes[genome].append(gene)
                annotations[genome][gene] = ['', '']
            elif (i % 4 == 2 or i % 4 == 3) and not empty:
                annotations[genome][gene][(i % 4) - 2] = row[0]

for cluster in clusters:
    annotationToLift = ['', '']
    for genome, gene in cluster:
        if annotations[genome][gene][0] != '' or annotations[genome][gene][1] != '':
            annotationToLift = annotations[genome][gene]
            break
    for genome, gene in cluster:
        annotations[genome][gene] = annotationToLift

maxNumGenes = len(max(genes.values(), key=len))
genomes = genes.keys()
for i, genome in enumerate(genomes):
    if i != 0:
        sys.stdout.write('\t')
    sys.stdout.write('"%s"' % genome)
sys.stdout.write("\n")

for i in xrange(maxNumGenes * 4):
    for j, genome in enumerate(genomes):
        if j != 0:
            sys.stdout.write('\t')
        if i >= len(genes[genome]) * 4 or i % 4 == 1:
            sys.stdout.write('""')
        else:
            gene = genes[genome][i / 4]
            if i % 4 == 0:
                sys.stdout.write('"%s"' % gene)
            elif i % 4 == 2 or i % 4 == 3:
                sys.stdout.write('"%s"' % annotations[genome][gene][i % 4 - 2])
    sys.stdout.write("\n")
