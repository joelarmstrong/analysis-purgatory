#!/usr/bin/env python
"""Testing the accuracy of clustering methods for splitting
homologous blocks into reference-orthologous blocks.

Requires:
- scikit-learn
- kmodes
- biopython
- sonLib
- pandas
"""
from argparse import ArgumentParser
from StringIO import StringIO
from collections import defaultdict
import itertools
import random
from sklearn.cluster import KMeans
from sklearn.preprocessing import OneHotEncoder
from kmodes.kmodes import KModes
from Bio import Phylo
from Bio.Phylo import TreeConstruction
from Bio.Phylo.BaseTree import Clade, Tree
import numpy as np
import pandas as pd
from sonLib.bioio import fastaRead
from simulator import BirthDeathSimulator, GeneralizedReversibleSimulator

cluster_methods = ['k-modes', 'k-means', 'neighbor-joining', 'upgma']
evaluation_methods = ['split-decomposition', 'none']

def seqs_to_columns(seqs, seq_order):
    """
    Transform a dict of sequences into a list of columns.

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
    """
    Take a nested list of DNA columns and convert them into a np feature array.
    """
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
    parser.add_argument('species_tree', help='species tree (newick format)')
    parser.add_argument('--duplication-rate',
                        type=float,
                        default=0.2,
                        help='Gene duplication rate')
    parser.add_argument('--loss-rate',
                        type=float,
                        default=0.05,
                        help='Gene loss rate')
    parser.add_argument('--num-columns',
                        type=int,
                        default=200,
                        help='Number of columns')
    parser.add_argument('--num-tests',
                        type=int,
                        default=100)
    return parser.parse_args()

def cluster_matrix(matrix, cluster_method):
    """
    Run a clustering method on a matrix and return the cluster assignments.
    """
    if cluster_method == 'k-means':
        return KMeans(n_clusters=2).fit_predict(matrix)
    elif cluster_method == 'k-modes':
        return KModes(n_clusters=2).fit_predict(matrix.todense())
    else:
        raise ArgumentError('Unknown cluster method: %s' % cluster_method)

def distance_matrix_from_columns(columns):
    """
    Get a distance matrix (as a np array) from a nested list of DNA columns.
    """
    num_seqs = len(columns[0])
    matrix = np.zeros([num_seqs, num_seqs], dtype=int)
    for column in columns:
        for i, entry_1 in enumerate(column):
            for j in xrange(i + 1, len(column)):
                entry_2 = column[j]
                if entry_1.lower() != entry_2.lower():
                    matrix[i, j] += 1
                    matrix[j, i] += 1
    return np.true_divide(matrix, len(columns))

def satisfies_four_point_criterion(matrix, split1, split2, relaxed=False,
                                   enforce_three_point=True):
    """Tests whether a split satisfies the d-split criterion of Bandelt
    and Dress 1992.

    The "relaxed" version is the same version that is in the paper,
    where the internal distance may be larger than one of the
    inter-split distances. Otherwise, it must be smaller than both.

    The "enforce_three_point" parameter determines whether to enforce
    the inter-vs-intra distance comparison even when one side of the
    split is a singleton. This isn't justified by the tree metric, but
    may work in practice.
    """
    if len(split1) < len(split2):
        # We call split1 the larger split, to simplify the
        # "enforce_three_point" logic
        split1, split2 = split2, split1

    for i, j in itertools.combinations(split1, 2):
        if len(split2) == 1 and enforce_three_point:
            intra = matrix[i, j]
            k = split2[0]
            inter = matrix[i, k] + matrix[j, k]
            if intra > inter:
                return False
        for k, l in itertools.combinations(split2, 2):
            intra = matrix[i, j] + matrix[k, l]
            inter1 = matrix[i, k] + matrix[j, l]
            inter2 = matrix[i, l] + matrix[j, k]
            if relaxed:
                if intra > inter1 and intra > inter2:
                    return False
            else:
                if intra > inter1 or intra > inter2:
                    return False

    return True

def is_good_split(cluster_assignments, columns, evaluation_method):
    """
    Run the "split-evaluation" method on the columns, and return True if it calls this split good.
    """
    assert all([i == 0 or i == 1 for i in cluster_assignments]), \
        "A valid split should only split into two partitions"
    if evaluation_method == 'none':
        return True
    elif evaluation_method == 'split-decomposition':
        distance_matrix = distance_matrix_from_columns(columns)
        split1 = [i for i, cluster in enumerate(cluster_assignments) if cluster == 0]
        split2 = [i for i, cluster in enumerate(cluster_assignments) if cluster == 1]
        return satisfies_four_point_criterion(distance_matrix, split1, split2)

def build_tree_top_down(columns, seq_names, cluster_method, evaluation_method):
    """
    Build a tree top-down using successive applications of a clustering method.
    """
    def is_finished(cluster, columns):
        return len(cluster) <= 2 or all([len(set(column)) == 1 for column in columns])
    def recurse(columns, seq_names):
        matrix = columns_to_matrix(columns, one_hot_encode=(cluster_method == 'k-means'))
        cluster_assignments = cluster_matrix(matrix, cluster_method)
        if not is_good_split(cluster_assignments, columns, evaluation_method):
            return Clade(clades=map(lambda x: Clade(name=x), seq_names))
        cluster0_indices = [i for i, cluster in enumerate(cluster_assignments) if cluster == 0]
        cluster1_indices = [i for i, cluster in enumerate(cluster_assignments) if cluster == 1]
        cluster0 = [seq_names[i] for i in cluster0_indices]
        cluster1 = [seq_names[i] for i in cluster1_indices]
        cluster0_columns = [[column[i] for i in cluster0_indices] for column in columns]
        if is_finished(cluster0, cluster0_columns):
            if len(cluster0) == 1:
                clade0 = Clade(name=cluster0[0])
            else:
                clade0 = Clade(clades=map(lambda x: Clade(name=x), cluster0))
        else:
            clade0 = Clade(clades=recurse(cluster0_columns, cluster0))

        cluster1_columns = [[column[i] for i in cluster1_indices] for column in columns]
        if is_finished(cluster1, cluster1_columns):
            if len(cluster1) == 1:
                clade1 = Clade(name=cluster1[0])
            else:
                clade1 = Clade(clades=map(lambda x: Clade(name=x), cluster1))
        else:
            clade1 = Clade(clades=recurse(cluster1_columns, cluster1))
        return (clade0, clade1)
    tree = Tree(Clade(clades=recurse(columns, seq_names)))
    return tree

def random_sequence(length):
    """
    Get a random DNA sequence of a certain length.
    """
    seq = []
    for _ in xrange(length):
        seq.append(random.choice(['A', 'a', 'C', 'c', 'G', 'g', 'T', 't']))
    return seq

def generate_gene_tree_and_sequences(gene_tree_sim, grt_sim, num_columns):
    """
    Get a random gene tree and the corresponding sequences of the leaves.

    The root sequence a random DNA sequence of length num_columns.
    """
    gene_tree = gene_tree_sim.generate()
    seqs = grt_sim.generate_leaf_sequences(gene_tree, random_sequence(num_columns))
    return gene_tree, seqs

def flatten_list(l):
    """
    Flatten a list.
    """
    return [item for sublist in l for item in sublist]

def build_tree_bottom_up(columns, seq_names, cluster_method, evaluation_method):
    """
    Build a tree using a NJ-esque clustering method.
    """
    distance_matrix = distance_matrix_from_columns(columns)
    # Biopython annoyingly (but understandably) wants the matrix in
    # lower triangular format, i.e. only everything below the diagonal
    triangular_matrix = [[entry for j,entry in enumerate(row) if j <= i] for i, row in enumerate(distance_matrix.tolist())]
    tree_constructor = TreeConstruction.DistanceTreeConstructor()
    distance_matrix = TreeConstruction._DistanceMatrix(seq_names, triangular_matrix)
    if cluster_method == 'neighbor-joining':
        tree = tree_constructor.nj(distance_matrix)
    elif cluster_method == 'upgma':
        tree = tree_constructor.upgma(distance_matrix)
    else:
        raise ArgumentError('Unrecognized bottom-up method: %s' % cluster_method)
    for internal_node in tree.get_nonterminals():
        split = [child for child in internal_node]
        if len(split) != 2:
            continue
        split1 = set([leaf.name for leaf in split[1].get_terminals()])
        leaf_names = [node.name for node in internal_node.get_terminals()]
        relevant_seq_names = [name for name in seq_names if name in leaf_names]
        relevant_indices = [i for i, name in enumerate(seq_names) if name in relevant_seq_names]
        relevant_columns = [[column[i] for i in relevant_indices] for column in columns]
        cluster_assignments = [int(seq_name in split1) for seq_name in relevant_seq_names]
        if not is_good_split(cluster_assignments, relevant_columns, evaluation_method):
            # Need to make this node into a multifurcation.
            internal_node.clades = flatten_list([[grandchild for grandchild in child] for child in split])

    return tree

def build_tree(seqs, cluster_method, evaluation_method, outgroups):
    """
    Build a tree using some clustering method and some split-evaluation method.
    """
    seq_names = seqs.keys()
    cols = seqs_to_columns(seqs, seq_names)
    if cluster_method in ['k-means', 'k-modes']:
        tree = build_tree_top_down(cols, seq_names, cluster_method, evaluation_method)
    else:
        tree = build_tree_bottom_up(cols, seq_names, cluster_method, evaluation_method)
    # workaround for biopython bug.
    for node in tree.find_clades():
        node.clades = list(node.clades)
    tree.root_with_outgroup(*[node for node in tree.get_terminals() if node.name in outgroups])
    return tree

def evaluate_tree(true_tree, test_tree):
    """
    Given a true tree and a test tree, give some stats on how well the test tree reflects the truth.
    """
    true_leaf_sets = set()
    true_splits = {}
    for internal_node in true_tree.get_nonterminals():
        leaf_set = frozenset([node.name for node in internal_node.get_terminals()])
        split_sets = frozenset([frozenset([node.name for node in node.get_terminals()]) for node in internal_node])
        true_leaf_sets.add(leaf_set)
        true_splits[leaf_set] = split_sets

    overcollapses = 0
    undercollapses = 0
    wrong_splits = 0
    mismatching_leaf_sets = 0
    perfect_splits = 0
    for internal_node in test_tree.get_nonterminals():
        leaf_set = frozenset([node.name for node in internal_node.get_terminals()])
        if leaf_set in true_leaf_sets:
            split_sets = frozenset([frozenset([node.name for node in node.get_terminals()]) for node in internal_node])
            if split_sets == true_splits[leaf_set]:
                perfect_splits += 1
            else:
                true_split = true_splits[leaf_set]
                if len(true_split) > len(split_sets):
                    overcollapses += 1
                elif len(true_split) < len(split_sets):
                    undercollapses += 1
                else:
                    wrong_splits += 1
        else:
            mismatching_leaf_sets += 1
    return { 'overcollapses': overcollapses,
             'undercollapses': undercollapses,
             'wrong_splits': wrong_splits,
             'mismatching_leaf_sets': mismatching_leaf_sets,
             'perfect_splits': perfect_splits }

def main():
    args = parse_args()
    species_tree = Phylo.read(StringIO(args.species_tree), 'newick')
    gene_tree_sim = BirthDeathSimulator(species_tree,
                                        args.duplication_rate,
                                        args.loss_rate)
    grt_sim = GeneralizedReversibleSimulator(0.25, 0.25, 0.25, 0.25,
                                             0.25, 0.25, 0.25, 0.25, 0.25)
    tree_evaluations = []
    for _ in xrange(args.num_tests):
        true_tree, leaf_seqs = generate_gene_tree_and_sequences(gene_tree_sim, grt_sim,
                                                                args.num_columns)
        for cluster_method in cluster_methods:
            for evaluation_method in evaluation_methods:
                # Choose the second child of the root as the outgroup for no good reason
                outgroups = [node.name for node in true_tree.root[1].get_terminals()]
                built_tree = build_tree(leaf_seqs, cluster_method, evaluation_method, outgroups)
                evaluation = evaluate_tree(true_tree, built_tree)
                evaluation['cluster_method'] = cluster_method
                evaluation['evaluation_method'] = evaluation_method
                tree_evaluations.append(evaluation)

    df = pd.DataFrame(tree_evaluations)
    print df.groupby(['cluster_method', 'evaluation_method']).sum()

if __name__ == '__main__':
    main()
