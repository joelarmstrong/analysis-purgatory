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
from sklearn import linear_model
from sklearn.feature_extraction import DictVectorizer
from sklearn.metrics import mean_squared_error

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
            datum = json.loads(match.group(2))
            # Now that we shrunk things, we only (probably) care about one feature
            datum = {k: datum[k] for k in datum if k in ['MaxMemory', 'maxFlowerSize']}
            # Give it extra headroom
            datum['MaxMemory'] = int(match.group(3)) * (1.0 + args.headroom)
            job = match.group(1)
            jobData[job].append(datum)

    for job, data in jobData.iteritems():
        vec = DictVectorizer()
        y = np.array([d.pop('MaxMemory') for d in data])
        X = vec.fit_transform(data)
        lm = linear_model.LinearRegression()
        lm.fit(X, y)
        # y without headroom adjustment
        true_y = y / (1.0 + args.headroom)
        reprediction = lm.predict(X)
        mse = mean_squared_error(true_y, reprediction)
        max_overestimate = max(reprediction - true_y)
        max_underestimate = max(true_y - reprediction)
        print 'Job %s: r^2: %s MSE: %s max overestimate: %sGiB max underestimate: %sGiB coef: %s intercept: %s' % \
            (job, lm.score(X, y), mse,
             max_overestimate/(1024*1024*1024), max_underestimate/(1024*1024*1024),
             dict(zip(vec.get_feature_names(), lm.coef_)), lm.intercept_)

if __name__ == '__main__':
    main()
