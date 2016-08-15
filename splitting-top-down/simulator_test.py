#!/usr/bin/env python
import unittest
import random
from hypothesis import given, assume
from hypothesis.strategies import text, builds, floats, sampled_from, composite, random_module, integers
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree, Clade
from StringIO import StringIO
from simulator import GeneralizedReversibleSimulator, BirthDeathSimulator

# random data strategies
probability = floats(min_value=0.0, max_value=1.0)
non_zero_probability = floats(min_value=0.01, max_value=1.0)
random_DNA = text(alphabet=['A', 'a', 'C', 'c', 'G', 'g', 'T', 't'])
# We need to use only somewhat realistic distances, because the
# matrix exponential is only approximate and becomes
# inaccurate at very high distances.
random_distance = floats(min_value=0.0, max_value=5)
@composite
def random_tree(draw, max_depth=5):
    root = draw(random_clade(max_depth=max_depth))
    root.branch_length = None
    return Tree(root)
@composite
def random_clade(draw, depth=0, max_depth=8):
    name = draw(text())
    branch_length = draw(random_distance)
    children = []
    if depth < max_depth:
        num_children = draw(integers(min_value=0, max_value=4))
        for _ in xrange(num_children):
            children.append(draw(random_clade(depth=depth+1, max_depth=max_depth)))
    return Clade(name=name, branch_length=branch_length, clades=children)

@composite
def randomGRT(draw):
    frac_a = draw(non_zero_probability)
    frac_c = draw(non_zero_probability)
    frac_g = draw(non_zero_probability)
    frac_t = draw(non_zero_probability)
    # Normalize the equilibrium frequencies
    sum_frac = frac_a + frac_c + frac_g + frac_t
    frac_a /= sum_frac
    frac_c /= sum_frac
    frac_g /= sum_frac
    frac_t /= sum_frac

    a_c = draw(non_zero_probability)
    a_g = draw(non_zero_probability)
    a_t = draw(non_zero_probability)
    c_g = draw(non_zero_probability)
    c_t = draw(non_zero_probability)
    g_t = draw(non_zero_probability)

    # Normalize the change parameters
    sum_change = 2*frac_a*frac_c*a_c + \
                 2*frac_a*frac_g*a_g + \
                 2*frac_a*frac_t*a_t + \
                 2*frac_c*frac_g*c_g + \
                 2*frac_c*frac_t*c_t + \
                 2*frac_g*frac_t*g_t

    a_c /= sum_change
    a_g /= sum_change
    a_t /= sum_change
    c_g /= sum_change
    c_t /= sum_change
    g_t /= sum_change

    sum_change = 2*frac_a*frac_c*a_c + \
                 2*frac_a*frac_g*a_g + \
                 2*frac_a*frac_t*a_t + \
                 2*frac_c*frac_g*c_g + \
                 2*frac_c*frac_t*c_t + \
                 2*frac_g*frac_t*g_t

    return GeneralizedReversibleSimulator(frac_a, frac_c, frac_g, a_c, a_g, a_t, c_g, c_t, g_t)

class GRTSimulatorTest(unittest.TestCase):
    char_to_frac_param = { 'A': 'frac_a', 'C': 'frac_c', 'G': 'frac_g' }
    chars_to_change_param = { ('A', 'C'): 'a_c', ('A', 'G'): 'a_g', ('A', 'T'): 'a_t',
                              ('C', 'G'): 'c_g', ('C', 'T'): 'c_t', ('G', 'T'): 'g_t' }
    default_params = { 'frac_a': 0.25, 'frac_c': 0.25, 'frac_g': 0.25,
                       'a_c': 0.25, 'a_g': 0.25, 'a_t': 0.25, 'c_g': 0.25, 'c_t': 0.25,
                       'g_t': 0.25 }
    default_sim = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25,
                                                 a_c=0.25, a_g=0.25, a_t=0.25, c_g=0.25, c_t=0.25,
                                                 g_t=0.25)

    def test_increasing_change_probability(self):
        """Increasing the change probability should be reflected in the probability."""
        initial_probability = self.default_sim.probability('A', 'T', 1.0)
        sim = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25,
                                             a_c=0.25, a_g=0.25, a_t=0.4, c_g=0.25, c_t=0.25, g_t=0.25)
        self.assertGreater(sim.probability('A', 'T', 1.0), initial_probability)

    def test_increasing_proportion(self):
        """Increasing the proportion of a character should be reflected in the probability."""
        for char in ['A', 'C', 'G']:
            params = self.default_params.copy()
            params[self.char_to_frac_param[char]] = 0.4
            other_chars = [c for c in ['A', 'C', 'G'] if c != char]
            for other_char in other_chars:
                params[self.char_to_frac_param[other_char]] = 0.2
            sim = GeneralizedReversibleSimulator(**params)
            initial_probability = self.default_sim.probability(char, char, 1.0)
            self.assertGreater(sim.probability(char, char, 1.0), initial_probability)

    @given(randomGRT(), sampled_from(['A', 'C', 'G', 'T']), random_distance)
    def test_probability_sums_to_1(self, sim, char, distance):
        """Test that the probability from a character to all characters sums to 1.0."""
        assume(distance > 0)
        total_probability = sim.probability(char, 'A', distance) + sim.probability(char, 'C', distance) + sim.probability(char, 'G', distance) + sim.probability(char, 'T', distance)
        self.assertAlmostEqual(total_probability, 1.0)

    @given(randomGRT(), random_DNA, random_distance)
    def test_mutate_gives_same_length_sequence(self, sim, sequence, distance):
        mutated = sim.mutate(sequence, distance)
        self.assertEqual(len(mutated), len(sequence))

    @given(randomGRT(), random_DNA)
    def test_generate_leaf_sequences_gives_same_length_sequence(self, sim, sequence):
        species_tree = Phylo.read(StringIO('((((HUMAN:0.006969, CHIMP:0.009727):0.025291, RHESUS:0.044568):0.11,(MOUSE:0.072818, RAT:0.081244):0.260342):0.023260,((DOG:0.07, CAT:0.07):0.087381,((PIG:0.06, COW:0.06):0.104728,HORSE:0.05):0.05):0.04);'), 'newick')
        leaf_sequences = sim.generate_leaf_sequences(species_tree, sequence)
        self.assertEqual(len(leaf_sequences), 10)
        self.assertTrue(all([leaf in leaf_sequences for leaf in ['HUMAN', 'CHIMP', 'RHESUS', 'MOUSE', 'RAT', 'DOG', 'CAT', 'PIG', 'COW', 'HORSE']]))
        self.assertTrue(all([len(leaf_sequence) == len(sequence) for leaf_sequence in leaf_sequences.values()]))

class BirthDeathSimulatorTest(unittest.TestCase):
    @given(random_tree(), probability, probability, random_module())
    def test_extinct_lineages_are_pruned(self, tree, duplication_rate, loss_rate, random_module):
        # Duplication rate should be reasonable (if it is high we get
        # an explosion in gene tree size)
        assume(duplication_rate < 0.2)
        seed = random.random()
        sim = BirthDeathSimulator(tree, duplication_rate, loss_rate)
        random.seed(seed)
        tree_with_extinctions = sim.generate(remove_extinct_lineages=False)
        random.seed(seed)
        tree_without_extinctions = sim.generate()
        names_with_extinctions = [node.name for node in tree_with_extinctions.get_terminals() if node != tree_with_extinctions.root]
        names_without_extinctions = [node.name for node in tree_without_extinctions.get_terminals() if node != tree_without_extinctions.root]
        self.assertEqual(names_without_extinctions,
                         [name for name in names_with_extinctions if 'extinct' not in name])

if __name__ == '__main__':
    unittest.main()
