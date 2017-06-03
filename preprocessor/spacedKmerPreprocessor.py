#!/usr/bin/env python
"""
Quick stab at replacing our very expensive preprocessor (basically
a self-alignment of all input genomes) with a much faster version that
just masks overrepresented seeds.

Just a prototype, not meant to be particularly fast.
"""
from sonLib.bioio import fastaRead
from argparse import ArgumentParser
from collections import Counter

def pack(word, pattern):
    """Return a packed word given a spaced seed pattern.

    >>> pack('actgac', [True, False, True, True, False, True])
    'atgc'
    """
    assert len(word) == len(pattern)
    ret = []
    for i, char in enumerate(word):
        if pattern[i]:
            ret.append(char)
    return "".join(ret)

def yieldSpacedKmers(seq, pattern):
    for i in xrange(len(seq) - len(pattern) + 1):
        yield pack(seq[i:i+len(pattern)], pattern)

def countSpacedKmers(seq, pattern):
    """Returns a Counter object containing the count of each seed.

    >>> c = countSpacedKmers('actgactg', [True, False, True])
    >>> c == Counter(['at', 'cg', 'ta', 'gc', 'at', 'cg'])
    True
    """
    return Counter(yieldSpacedKmers(seq, pattern))

def parse_pattern(pattern):
    """Convert a string of 1s and 0s into a pattern of Trues and Falses.

    >>> parse_pattern('1010010')
    [True, False, True, False, False, True, False]
    """
    return map(lambda x: True if x == '1' else False, pattern)

def parse_args():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('fasta')
    parser.add_argument('--pattern', type=parse_pattern,
                        default='1110100110010101111')
    return parser.parse_args()

def main():
    opts = parse_args()
    count = Counter()
    for header, seq in fastaRead(opts.fasta):
        print 'read'
        seq = seq.lower()
        print 'lowered'
        seq = seq.translate(None, 'n')
        print 'filtered'
        kmers = countSpacedKmers(seq, opts.pattern)
        print header, kmers.most_common(10)
        count += kmers
    print count.most_common(10)
    with open('dump', 'w') as f:
        for key, val in count.items():
            f.write('%s\t%s\n' % (key, val))

if __name__ == '__main__':
    main()
