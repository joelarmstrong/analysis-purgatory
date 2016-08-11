#!/usr/bin/env python
from hypothesis import given, assume
from hypothesis.strategies import builds, floats, sampled_from, composite
import unittest
from simulator import GeneralizedReversibleSimulator

# random data strategies
probability = floats(min_value=0.0, max_value=1.0)
non_zero_probability = floats(min_value=0.01, max_value=1.0)
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

    @given(randomGRT(),
           sampled_from(['A', 'C', 'G', 'T']),
           # We need to use only somewhat realistic distances, because the
           # matrix exponential is only approximate and becomes
           # inaccurate at very high distances.
           floats(min_value=0.0, max_value=500000))
    def test_probability_sums_to_1(self, sim, char, distance):
        """Test that the probability from a character to all characters sums to 1.0."""
        assume(distance > 0)
        total_probability = sim.probability(char, 'A', distance) + sim.probability(char, 'C', distance) + sim.probability(char, 'G', distance) + sim.probability(char, 'T', distance)
        self.assertAlmostEqual(total_probability, 1.0)

if __name__ == '__main__':
    unittest.main()
