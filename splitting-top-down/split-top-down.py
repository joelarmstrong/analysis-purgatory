#!/usr/bin/env python
"""
Testing the accuracy of a top-down clustering method using k-means
clustering, k = 2, for splitting homologous blocks into
reference-orthologous blocks.
"""
from argparse import ArgumentParser
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from kmodes.kmodes import KModes
import numpy as np
from sonLib.bioio import fastaRead

def seqs_to_columns(seqs, seq_order):
    """Transform a dict of sequences into a list of columns.

    Each column is represented by a list of entries. The order of each
    sequence within the entries is the same as the order in the
    parameter seq_order.
    """
    assert len(seq_order) == len(seqs.keys()), \
        "'seq_order' has more or fewer entries than 'seqs'"
    assert all([seq_name in seqs for seq_name in seq_order]), \
        "'seq_order' refers to a sequence not present in 'seqs'"
    assert len(set([len(seq) for seq in seqs.values()])), \
        "All sequences must have the same length"
    columns = []
    for seq_name in seq_order:
        seq = seqs[seq_name]
        if len(columns) == 0:
            columns = [[] for _ in xrange(len(seq))]
        for i, char in enumerate(seq):
            columns[i].append(char)
    return columns

def columns_to_matrix(cols, one_hot_encode=True):
    def nuc_to_number(nucleotide):
        nucleotide = nucleotide.lower()
        if nucleotide == 'a':
            return 0
        elif nucleotide == 'c':
            return 1
        elif nucleotide == 'g':
            return 2
        elif nucleotide == 't':
            return 3
        else:
            return 4
    transformed_cols = [map(nuc_to_number, col) for col in cols]
    raw_matrix = np.matrix(transformed_cols).transpose()
    encoder = OneHotEncoder()
    encoded_matrix = encoder.fit_transform(raw_matrix)
    return encoded_matrix

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('fasta', help='fasta file (all sequences must have the '
                        'same length')
    parser.add_argument('species-tree', help='species tree')
    parser.add_argument('reference', help='reference species')
    parser.add_argument('--cluster-method',
                        choices=['k-means', 'k-modes'],
                        default='k-means',
                        help='Clustering method to use')
    return parser.parse_args()

def cluster_matrix(matrix, cluster_method):
    if cluster_method == 'k-means':
        return KMeans(n_clusters=2).fit_predict(matrix)
    elif cluster_method == 'k-modes':
        return KModes(n_clusters=2).fit_predict(matrix.todense())
    else:
        raise ArgumentError('Unknown cluster method: %s' % cluster_method)

def main():
    args = parse_args()
    seqs = dict(fastaRead(args.fasta))
    seq_names = seqs.keys()
    print repr(seq_names)
    cols = seqs_to_columns(seqs, seq_names)
    matrix = columns_to_matrix(cols, one_hot_encode=(args.cluster_method == 'k-means'))
    print cluster_matrix(matrix, args.cluster_method)

if __name__ == '__main__':
    main()
