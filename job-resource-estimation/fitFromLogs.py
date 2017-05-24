#!/usr/bin/env python
"""
Generates a polynomial fit to the memory requirements of jobs
given a set of toil logs.

By default, it generates a quadratic fit.
"""
import re
import json
from argparse import ArgumentParser
from collections import defaultdict

import numpy as np
import cvxpy
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import mean_squared_error

def fit_lm_never_underestimate(X, y):
    """Fit a linear model to the data, subject to the constraint that the
    linear model never underestimates the training points."""
    n = X.shape[1]    # Number of features
    m = X.shape[0]    # Number of examples
    B = cvxpy.Variable(n)
    # Must always be above the training points
    constraints = [(X*B)[i] >= y[i] for i in xrange(m)]
    # Minimize absolute error
    objective = cvxpy.Minimize(cvxpy.sum_entries(X*B - y))
    problem = cvxpy.Problem(objective, constraints)
    problem.solve()
    return B.value

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('logFile')
    parser.add_argument('--headroom', help='Add an extra x ratio of headroom',
                        type=float, default=0.1)
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    jobData = defaultdict(list)
    with open(args.logFile) as f:
        for line in f:
            match = re.search(r'Max memory used for job (.+) on JSON features (.+): ([0-9]+)', line)
            assert match is not None
            job = match.group(1)
            datum = json.loads(match.group(2))
            # Now that we shrunk things, we only (probably) care about one feature
            datum = {k: datum[k] for k in datum if k in ['MaxMemory', 'maxFlowerSize']}
            # Give it extra headroom
            datum['MaxMemory'] = int(match.group(3)) * (1.0 + args.headroom)
            jobData[job].append(datum)

    for job, data in jobData.iteritems():
        vec = DictVectorizer(sparse=False)
        y = np.array([d.pop('MaxMemory') for d in data])
        X = vec.fit_transform(data)
        # Extend the data so we can fit an intercept
        X = np.concatenate([np.ones((X.shape[0], 1)), X], axis=1)
        B = fit_lm_never_underestimate(X, y)
        # y without headroom adjustment
        reprediction = X.dot(B)
        true_y = (y / (1.0 + args.headroom)).reshape(reprediction.shape)
        mse = mean_squared_error(true_y, reprediction)
        max_overestimate = max(reprediction - true_y)
        max_underestimate = max(true_y - reprediction)
        print 'Job %s: MSE: %s max overestimate: %sGiB max underestimate: %sGiB coef: %s' % \
            (job, mse,
             max_overestimate/(1024*1024*1024), max_underestimate/(1024*1024*1024),
             dict(zip(['intercept'] + vec.get_feature_names(), B)))

if __name__ == '__main__':
    main()
