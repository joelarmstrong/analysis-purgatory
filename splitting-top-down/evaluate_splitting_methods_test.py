#!/usr/bin/env python
import unittest
import numpy as np
from evaluate_splitting_methods import get_d_splits

class DSplitTest(unittest.TestCase):
    def setUp(self):
        self.distance_matrix = np.array([[ 0.,  4.,  5.,  7., 13.,  8.,  6.],
                                         [ 4.,  0.,  1.,  3.,  9., 12., 10.],
                                         [ 5.,  1.,  0.,  2.,  8., 13., 11.],
                                         [ 7.,  3.,  2.,  0.,  6., 11., 13.],
                                         [13.,  9.,  8.,  6.,  0.,  5.,  7.],
                                         [ 8., 12., 13., 11.,  5.,  0.,  2.],
                                         [ 6., 10., 11., 13.,  7.,  2.,  0.]])

    def test_d_split(self):
        d_splits = get_d_splits(self.distance_matrix, relaxed=True)
        self.assertEqual(d_splits, set([(frozenset([0, 1, 2, 3]), frozenset([4, 5, 6])),
                                        (frozenset([0, 5, 6]), frozenset([1, 2, 3, 4])),
                                        (frozenset([0, 1, 2, 6]), frozenset([3, 4, 5])),
                                        (frozenset([0, 1, 5, 6]), frozenset([2, 3, 4]))]))

if __name__ == '__main__':
    unittest.main()
