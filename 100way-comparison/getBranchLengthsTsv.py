#!/usr/bin/env python
from argparse import ArgumentParser
from sonLib.bioio import popenCatch
parser = ArgumentParser()
parser.add_argument('halPath')
parser.add_argument('outputTsv')

opts = parser.parse_args()

genomes = popenCatch("halStats --genomes %s" % (opts.halPath)).split()
outFile = open(opts.outputTsv, 'w')
for genome in genomes:
    branchLength = float(popenCatch("halStats --branchLength %s"))
    outFile.write("%s\t%f\n" % (genome, branchLength))
