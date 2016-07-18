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
        clusters.append(tuple(cluster))

# map from genome name -> gene names (ordered)
genes = defaultdict(list)
# map from genome name -> gene name -> annotation (which may be empty)
annotations = defaultdict(dict)
# map from genome name -> gene name -> url
urls = defaultdict(dict)
# annotations per cluster
clusterAnnotations = defaultdict(list)

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
            elif i % 4 == 1:
                urls[genome][gene] = row[0]
            elif (i % 4 == 2 or i % 4 == 3) and not empty:
                annotations[genome][gene][(i % 4) - 2] = row[0]

for cluster in clusters:
    clusterAnnotations[cluster] = []
    for genome, gene in cluster:
        if annotations[genome][gene][0] != '' or annotations[genome][gene][1] != '':
            clusterAnnotations[cluster].append(annotations[genome][gene])

maxNumGenes = len(max(genes.values(), key=len))
genomes = genes.keys()
# write header
sys.stdout.write('"Description"')
for genome in genomes:
    sys.stdout.write('\t"%s ID"\t"%s URL"' % (genome, genome))
sys.stdout.write("\n")

for cluster in clusters:
    sys.stdout.write('"%s"' % (','.join(['|'.join(x) for x in clusterAnnotations[cluster]])))
    for genome in genomes:
        clusterRecord = [x for x in cluster if x[0] == genome]
        if len(clusterRecord) >= 1:
            # Multiple "genes" (transcripts) may be in the same genome
            # in the same cluster. Show the URL for only the first
            # gene but display both transcript IDs.
            firstGene = clusterRecord[0][1]
            gene = ','.join([record[1] for record in clusterRecord])
            sys.stdout.write('\t"%s"\t"%s"' % (gene, urls[genome][firstGene]))
        else:
            sys.stdout.write('\t""\t""')
    sys.stdout.write("\n")
