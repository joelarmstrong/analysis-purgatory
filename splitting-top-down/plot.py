#!/usr/bin/env python
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    sns.set_style('ticks')

    df = pd.read_csv('results.csv')
    df['fraction_good_splits'] = (df['perfect_splits']/df['fraction_perfect_splits'] - df['mismatching_leaf_sets'] - df['flipped_splits']) / (df['perfect_splits']/df['fraction_perfect_splits'])
    grid = sns.FacetGrid(df, size=5, row='evaluation_method',
                         row_order=['none', 'relaxed-split-decomposition', 'split-decomposition'],
                         hue_order=['k-means', 'k-modes', 'maximum-likelihood', 'upgma', 'neighbor-joining', 'guided-neighbor-joining', 'split-decomposition'],
                         aspect=16.0/9, legend_out=True)
    grid.map(sns.boxplot, 'loss_rate', 'fraction_good_splits', 'cluster_method', palette='colorblind', hue_order=['k-means', 'k-modes', 'maximum-likelihood', 'upgma', 'neighbor-joining', 'guided-neighbor-joining', 'split-decomposition']).set_axis_labels('Loss rate (as fraction of substitution rate)', 'Fraction of true splits correctly split or unresolved')
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 1.5))
    sns.plt.savefig('varying_loss_rate.pdf', bbox_extra_artists=(legend,), bbox_inches='tight')

    grid = sns.FacetGrid(df, size=5, row='evaluation_method',
                         row_order=['none', 'relaxed-split-decomposition', 'split-decomposition'],
                         hue_order=['k-means', 'k-modes', 'maximum-likelihood', 'upgma', 'neighbor-joining', 'guided-neighbor-joining', 'split-decomposition'],
                         aspect=16.0/9, legend_out=True)
    grid.map(sns.boxplot, 'duplication_rate', 'fraction_good_splits', 'cluster_method', palette='colorblind', hue_order=['k-means', 'k-modes', 'maximum-likelihood', 'upgma', 'neighbor-joining', 'guided-neighbor-joining', 'split-decomposition']).set_axis_labels('Duplication rate (as fraction of substitution rate)', 'Fraction of true splits correctly split or unresolved')
    legend = plt.legend(loc='center left', bbox_to_anchor=(1, 1.5))
    sns.plt.savefig('varying_duplication_rate.pdf', bbox_extra_artists=(legend,), bbox_inches='tight')


    print df.groupby(['cluster_method', 'evaluation_method']).sum().to_csv()

if __name__ == '__main__':
    main()
