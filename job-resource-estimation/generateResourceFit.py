#!/usr/bin/env python
from argparse import ArgumentParser

import bs4
from bs4 import BeautifulSoup
import pandas as pd
import numpy as np

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('statsXMLs', nargs='+')
    parser.add_argument('--lengths', required=True, type=int, nargs='+')
    parser.add_argument('--names', required=True, nargs='+')
    args = parser.parse_args()
    if len(args.names) != len(args.statsXMLs) \
       or len(args.lengths) != len(args.statsXMLs):
        raise RuntimeError("--names, --lengths, and statsXMLs must have the "
                           "same number of arguments")
    return args

def main():
    args = parse_args()
    data = []
    for name, length, statsXML in zip(args.names, args.lengths, args.statsXMLs):
        datum = { 'Name': name,
                  'TotalLength': length }
        soup = BeautifulSoup(open(statsXML), "xml")
        for i in soup.target_types.children:
            if isinstance(i, bs4.element.Tag):
                datum[i.name] = float(i['max_memory']) * 1024
        data.append(datum)
    df = pd.DataFrame(data)
    for col in df:
        if col not in ['Name', 'TotalLength']:
            print col
            print np.polyfit(df['TotalLength'], df[col], 2)

if __name__ == '__main__':
    main()
