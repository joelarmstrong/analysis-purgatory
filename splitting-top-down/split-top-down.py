#!/usr/bin/env python
"""
Testing the accuracy of a top-down clustering method using k-means
clustering, k = 2, for splitting homologous blocks into
reference-orthologous blocks.

Requires:
- scikit-learn
- kmodes
- biopython
- sonLib
"""
from argparse import ArgumentParser
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from kmodes.kmodes import KModes
from Bio import Phylo
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
    parser.add_argument('species_tree', help='species tree (newick format)')
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

def build_tree_topdown(columns, seq_names, cluster_method):
    def is_finished(cluster, columns):
        return len(cluster) <= 2 or all([len(set(column)) == 1 for column in columns])
    def recurse(columns, seq_names):
        matrix = columns_to_matrix(columns, one_hot_encode=(cluster_method == 'k-means'))
        cluster_assignments = cluster_matrix(matrix, cluster_method)
        cluster0_indices = [i for i, cluster in enumerate(cluster_assignments) if cluster == 0]
        cluster1_indices = [i for i, cluster in enumerate(cluster_assignments) if cluster == 1]
        cluster0 = [seq_names[i] for i in cluster0_indices]
        cluster1 = [seq_names[i] for i in cluster1_indices]
        cluster0_columns = [[column[i] for i in cluster0_indices] for column in columns]
        if is_finished(cluster0, cluster0_columns):
            clade0 = Phylo.BaseTree.Clade(clades=map(lambda x: Phylo.BaseTree.Clade(name=x), cluster0))
        else:
            clade0 = Phylo.BaseTree.Clade(clades=recurse(cluster0_columns, cluster0))

        cluster1_columns = [[column[i] for i in cluster1_indices] for column in columns]
        if is_finished(cluster1, cluster1_columns):
            clade1 = Phylo.BaseTree.Clade(clades=map(lambda x: Phylo.BaseTree.Clade(name=x), cluster1))
        else:
            clade1 = Phylo.BaseTree.Clade(clades=recurse(cluster1_columns, cluster1))
        return (clade0, clade1)
    tree = Phylo.BaseTree.Tree(Phylo.BaseTree.Clade(clades=recurse(columns, seq_names)))
    return tree

def main():
    args = parse_args()
    species_tree = Phylo.parse(args.species_tree, 'newick')
    seqs = dict(fastaRead(args.fasta))
    seq_names = seqs.keys()
    cols = seqs_to_columns(seqs, seq_names)
    tree = build_tree_topdown(cols, seq_names, args.cluster_method)
    Phylo.draw_ascii(tree)

if __name__ == '__main__':
    main()
