#!/usr/bin/env python
from scipy import linalg
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree
from StringIO import StringIO
import random
import sys
import numpy as np

def get_parent(tree, child_clade):
    """
    Get the parent of a node in a Bio.Phylo Tree.
    """
    if child_clade == tree.root:
        return None
    node_path = tree.root.get_path(child_clade)
    if len(node_path) < 2:
        return tree.root
    return node_path[-2]

def prune_lineage(tree, lineage):
    last_valid_lineage = None
    while lineage != None:
        parent = get_parent(tree, lineage)
        if parent is None:
            if len(lineage.clades) == 0:
                if last_valid_lineage is not None and last_valid_lineage.is_terminal():
                    return Tree(Clade(clades=[last_valid_lineage]))
                else:
                    return Tree(last_valid_lineage)
            else:
                assert len(lineage.clades) == 1
                if lineage.clades[0].is_terminal():
                    return Tree(Clade(clades=[lineage.clades[0]]))
                else:
                    return Tree(lineage.clades[0])
        assert lineage in parent
        parent.clades.remove(lineage)
        if len(parent.clades) >= 1:
            if len(parent.clades) == 1 and last_valid_lineage is None:
                last_valid_lineage = parent.clades[0]
            elif last_valid_lineage is not None:
                parent.clades.append(last_valid_lineage)
                return tree
        if len(parent.clades) <= 1:
            lineage = parent
        else:
            break
    return tree

def prune_lineages(tree, lineages_to_prune):
    """
    Prune a set of lineages from a Bio.Phylo Tree.
    """
    for lineage in lineages_to_prune:
        tree = prune_lineage(tree, lineage)
    return tree

class BirthDeathSimulator:
    """A simulator for generating gene trees using a birth-death process."""
    def __init__(self, species_tree, duplication_rate, extinction_rate):
        self.species_tree = species_tree
        self.duplication_rate = duplication_rate
        self.extinction_rate = extinction_rate

    def generate(self, dt=0.001, remove_extinct_lineages=True):
        self.max_lineage_num = 0
        def recurse(child_species, remaining_length, lineage_num, extinct_lineages):
            for time in np.arange(0, remaining_length, dt):
                # Convert from np.float64 to native python float,
                # because Bio.Phylo will flip out if float64's are
                # used as branch lengths
                time = float(time)
                if time + dt > remaining_length:
                    # Make sure the last step is for whatever time is
                    # still remaining, rather than just the full dt
                    timestep = remaining_length - time
                else:
                    timestep = dt
                duplication_occurred = random.random() < self.duplication_rate * timestep
                extinction_occurred = random.random() < self.extinction_rate * timestep
                # If both a duplication and an extinction
                # occurred, we assume only the extinction
                # happened, for simplicity.
                if extinction_occurred:
                    # Lineage is gone
                    extinction = Clade(name="extinct.%s" % lineage_num, branch_length=time)
                    extinct_lineages.append(extinction)
                    return extinction
                elif duplication_occurred:
                    # Lineage splits, add new lineage
                    self.max_lineage_num += 1
                    new_lineage_num = self.max_lineage_num
                    return Clade(branch_length=time,
                                 clades=[recurse(child_species, remaining_length - time, lineage_num, extinct_lineages),
                                         recurse(child_species, remaining_length - time, new_lineage_num, extinct_lineages)])
            if child_species.is_terminal():
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (child_species.name, lineage_num))
            else:
                return Clade(branch_length=remaining_length,
                             name="%s.%s" % (child_species.name, lineage_num),
                             clades=[recurse(child, child.branch_length, lineage_num, extinct_lineages) for child in child_species])
        extinct_lineages = []
        gene_tree = Tree(Clade(name='root', clades=[recurse(child, child.branch_length, 0, extinct_lineages) for child in self.species_tree.root]))
        if remove_extinct_lineages:
            gene_tree = prune_lineages(gene_tree, extinct_lineages)
        return gene_tree

class GeneralizedReversibleSimulator:
    """
    Simulate the evolution of DNA according to a generalized time-reversible process.
    """
    def __init__(self, frac_a, frac_c, frac_g, a_c, a_g, a_t, c_g, c_t, g_t):
        frac_t = 1.0 - frac_a - frac_c - frac_g
        self.rate_matrix = np.array([[0.0,          frac_c * a_c, frac_g * a_g, frac_t * a_t],
                                     [frac_a * a_c, 0.0,          frac_g * c_g, frac_t * c_t],
                                     [frac_a * a_g, frac_c * c_g, 0.0,          frac_t * g_t],
                                     [frac_a * a_t, frac_c * c_t, frac_g * g_t, 0.0         ]])
        # Calculate the diagonals
        row_sums = self.rate_matrix.sum(axis=1)
        self.rate_matrix[0, 0] = -row_sums[0]
        self.rate_matrix[1, 1] = -row_sums[1]
        self.rate_matrix[2, 2] = -row_sums[2]
        self.rate_matrix[3, 3] = -row_sums[3]

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

    def generate_leaf_sequences(self, tree, starting_sequence):
        root = tree.root
        node_to_sequence = { root: starting_sequence }
        def recurse(node, sequence, node_to_sequence):
            assert node_to_sequence[node] == sequence
            for child in node:
                child_sequence = self.mutate(sequence, child.branch_length)
                node_to_sequence[child] = child_sequence
                recurse(child, child_sequence, node_to_sequence)
        recurse(root, starting_sequence, node_to_sequence)
        return dict((node.name, sequence) for (node, sequence) in node_to_sequence.iteritems() if node.is_terminal())

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
            raise RuntimeError("Character not in {A,C,G,T}")

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
    mutator = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25,
                                             a_c=0.25, a_g=0.25, a_t=0.25, c_g=0.25, c_t=0.25, g_t=0.25)
    print mutator.mutate('ACGTACGTACGT', 1.0)
