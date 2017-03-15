#!/usr/bin/env python
"""
Generates a polynomial fit to the memory requirements of jobs
given a set of jobTreeStats XMLs.

By default, it generates a quadratic fit.
"""
from argparse import ArgumentParser

import bs4
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('statsXMLs', nargs='+')
    parser.add_argument('--lengths', required=True, type=int, nargs='+')
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--headroom', help='Add an extra x ratio of headroom',
                        type=float, default=0.1)
    args = parser.parse_args()
    if len(args.lengths) != len(args.statsXMLs):
        raise RuntimeError("--lengths and statsXMLs must have the "
                           "same number of arguments")
    return args

def evaluate(poly, x):
    """Evaluates a polynomial (provided as a pd.Series) at x."""
    X = [x**degree for degree in reversed(xrange(len(poly)))]
    return sum(poly*X)

def main():
    args = parse_args()
    data = []
    for length, statsXML in zip(args.lengths, args.statsXMLs):
        datum = { 'TotalLength': length }
        soup = BeautifulSoup(open(statsXML), "xml")
        for i in soup.target_types.children:
            if isinstance(i, bs4.element.Tag):
                datum[i.name] = float(i['max_memory']) * 1024
                # Give it extra headroom
                datum[i.name] += datum[i.name] * args.headroom
        data.append(datum)
    df = pd.DataFrame(data)
    for col in df:
        if col not in ['TotalLength']:
            df2 = df[['TotalLength', col]].dropna(axis=0)
            poly = np.polyfit(df2['TotalLength'], df2[col], args.degree)
            print poly
            print '%s: simMammals: %sG, nematodes: %sG, euarchontoglires: %sG' % \
                (col, evaluate(poly, 2400000)/(1024*1024*1024),
                 evaluate(poly, 429000000)/(1024*1024*1024),
                 evaluate(poly, 8279009913)/(1024*1024*1024))

if __name__ == '__main__':
    main()
