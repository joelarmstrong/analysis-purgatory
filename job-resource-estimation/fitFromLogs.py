#!/usr/bin/env python
"""
Generates a polynomial fit to the memory requirements of jobs
given a set of toil logs.

By default, it generates a quadratic fit.
"""
import re
from argparse import ArgumentParser

import pandas as pd
import numpy as np
import seaborn as sns

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('logFile')
    parser.add_argument('--degree', type=int, default=2)
    parser.add_argument('--headroom', help='Add an extra x ratio of headroom',
                        type=float, default=0.1)
    args = parser.parse_args()
    return args

def evaluate(poly, x):
    """Evaluates a polynomial (provided as a pd.Series) at x."""
    X = [x**degree for degree in reversed(xrange(len(poly)))]
    return sum(poly*X)

def main():
    args = parse_args()
    data = []
    with open(args.logFile) as f:
        for line in f:
            match = re.search(r'Max memory used for job (.+) on input size ([0-9]+): ([0-9]+)', line)
            assert match is not None
            datum = { 'Length': int(match.group(2)),
                      match.group(1): int(match.group(3)) }
            # Give it extra headroom
            datum[match.group(1)] += datum[match.group(1)] * args.headroom
            data.append(datum)

    df = pd.DataFrame(data)
    for col in df:
        if col not in ['Length']:
            df2 = df[['Length', col]].dropna(axis=0)
            poly = np.polyfit(df2['Length'], df2[col], args.degree)
            print poly
            print col
            sns.regplot(data=df2, x='Length', y=col, order=args.degree, ci=None)
            predicted = np.polyval(poly, df2['Length'])
            residuals = predicted - df2[col]
            # Plot residuals
            df3 = pd.DataFrame({'predicted':predicted, 'residuals': residuals})
            sns.jointplot(data=df3, x='predicted', y='residuals')
            sns.plt.show()

if __name__ == '__main__':
    main()
