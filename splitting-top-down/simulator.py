#!/usr/bin/env python
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from StringIO import StringIO
import random
import sys
import numpy as np

class BirthDeathSimulator:
    """A simulator for generating gene trees using a birth-death process."""
    def __init__(self, species_tree, duplication_rate, extinction_rate):
        self.species_tree = species_tree
        self.duplication_rate = duplication_rate
        self.extinction_rate = extinction_rate

    def generate(self, dt=0.001):
        self.max_lineage_num = 0
        def recurse(root_species, remaining_length, lineage_num):
            for time in np.arange(0, remaining_length, dt):
                # FIXME: will do the full dt on last iteration
                duplication_occurred = random.random() < self.duplication_rate * dt
                extinction_occurred = random.random() < self.extinction_rate * dt
                # If both a duplication and an extinction
                # occurred, we assume only the extinction
                # happened, for simplicity.
                if extinction_occurred:
                    # Lineage is gone
                    return Clade(name="extinct.%s" % lineage_num, branch_length=time)
                elif duplication_occurred:
                    # Lineage splits, add new lineage
                    self.max_lineage_num += 1
                    new_lineage_num = self.max_lineage_num
                    return Clade(branch_length=time,
                                 clades=[recurse(root_species, remaining_length - time, lineage_num),
                                         recurse(root_species, remaining_length - time, new_lineage_num)])
            if root_species.is_terminal():
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (root_species.name, lineage_num))
            else:
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (root_species.name, lineage_num),
                             clades=[recurse(child, child.branch_length, lineage_num) for child in root_species])
        return Tree(Clade(name='root', clades=[recurse(child, child.branch_length, 0) for child in self.species_tree.root]))

if __name__ == '__main__':
    species_tree = Phylo.read(StringIO(sys.argv[1]), 'newick')
    duplication_rate = 0.1
    extinction_rate = 0.05
    sim = BirthDeathSimulator(species_tree, duplication_rate, extinction_rate)
    Phylo.draw_ascii(sim.generate())
