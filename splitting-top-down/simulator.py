#!/usr/bin/env python
from scipy import linalg
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
        def recurse(child_species, remaining_length, lineage_num):
            for time in np.arange(0, remaining_length, dt):
                # Convert from np.float64 to native python float,
                # because Bio.Phylo will flip out if float64's are
                # used as branch lengths
                time = float(time)
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
                                 clades=[recurse(child_species, remaining_length - time, lineage_num),
                                         recurse(child_species, remaining_length - time, new_lineage_num)])
            if child_species.is_terminal():
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (child_species.name, lineage_num))
            else:
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (child_species.name, lineage_num),
                             clades=[recurse(child, child.branch_length, lineage_num) for child in child_species])
        return Tree(Clade(name='root', clades=[recurse(child, child.branch_length, 0) for child in self.species_tree.root]))

class GeneralizedReversibleSimulator:
    def __init__(self, frac_a, frac_c, frac_g, frac_t, a_c, a_g, a_t, c_g, c_t, g_t):
        self.rate_matrix = np.array([[0.0,          frac_c * a_c, frac_g * a_g, frac_t * a_t],
                                     [frac_a * a_c, 0.0,          frac_g * c_g, frac_t * c_t],
                                     [frac_a * a_g, frac_c * c_g, 0.0,          frac_t * g_t],
                                     [frac_a * a_t, frac_c * c_t, frac_g * g_t, 0.0         ]])
        # Calculate the diagonals
        row_sums = self.rate_matrix.sum(axis=1)
        self.rate_matrix[0, 0] = 1.0 - row_sums[0]
        self.rate_matrix[1, 1] = 1.0 - row_sums[1]
        self.rate_matrix[2, 2] = 1.0 - row_sums[2]
        self.rate_matrix[3, 3] = 1.0 - row_sums[3]

    def mutate(self, seq, distance):
        matrix = self.compute_probability_matrix(distance)
        matrix_rowsums = matrix.cumsum(axis=1)
        new_seq = []
        for char in seq:
            rand = random.random()
            new_char = None
            for i, cumsum in enumerate(matrix_rowsums[self.char_to_index(char)]):
                if rand <= cumsum:
                    new_char = self.index_to_char(i)
                    break
            assert new_char is not None
            new_seq.append(new_char)
        return "".join(new_seq)

    def char_to_index(self, char):
        char = char.lower()
        if char == 'a':
            return 0
        elif char == 'c':
            return 1
        elif char == 'g':
            return 2
        elif char == 't':
            return 3
        else:
            raise ArgumentError("Character not in {A,C,G,T}")

    def index_to_char(self, index):
        assert index >= 0 and index < 4
        return ['A', 'C', 'G', 'T'][index]

    def probability(self, from_char, to_char, distance, precomputed_matrix=None):
        if precomputed_matrix is not None:
            matrix = precomputed_matrix
        else:
            matrix = self.compute_probability_matrix(distance)
        from_index = self.char_to_index(from_char)
        to_index = self.char_to_index(to_char)
        return matrix[from_index, to_index]

    def compute_probability_matrix(self, distance):
        return linalg.expm(self.rate_matrix * distance)

if __name__ == '__main__':
    species_tree = Phylo.read(StringIO(sys.argv[1]), 'newick')
    duplication_rate = 0.1
    extinction_rate = 0.05
    print species_tree
    sim = BirthDeathSimulator(species_tree, duplication_rate, extinction_rate)
    tree = sim.generate()
    print tree
    Phylo.draw_ascii(tree)
    mutator = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25, frac_t=0.25,
                                             a_c=0.25, a_g=0.25, a_t=0.25, c_g=0.25, c_t=0.25, g_t=0.25)
    print mutator.mutate('ACGTACGTACGT', 1.0)
