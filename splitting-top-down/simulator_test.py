#!/usr/bin/env python
import unittest
from simulator import GeneralizedReversibleSimulator

class grttest(unittest.TestCase):
    def test_waffles(self):
        sim = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25, frac_t=0.25,
                                             a_c=0.25, a_g=0.25, a_t=0.25, c_g=0.25, c_t=0.25, g_t=0.25)
        initial_probability = sim.probability('A', 'T', 1.0)
        sim2 = GeneralizedReversibleSimulator(frac_a=0.25, frac_c=0.25, frac_g=0.25, frac_t=0.25,
                                             a_c=0.25, a_g=0.25, a_t=0.4, c_g=0.25, c_t=0.25, g_t=0.25)
        self.assertTrue(sim2.probability('A', 'T', 1.0)>initial_probability
)

        
if __name__ == '__main__':
    unittest.main()
