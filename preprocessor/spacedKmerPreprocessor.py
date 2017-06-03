#!/usr/bin/env python
"""
Quick stab at replacing our very expensive preprocessor (basically
a self-alignment of all input genomes) with a much faster version that
just masks overrepresented seeds.

Just a prototype, not meant to be particularly fast.
"""
from sonLib.bioio import fastaRead
from argparse import ArgumentParser
from collections import defaultdict

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
    """Yield spaced k-mers (and whether they are masked) according to the match/don't care pattern.

    k-mers containing Ns are skipped, and are considered masked if any bases are masked."""
    for i in xrange(len(seq) - len(pattern) + 1):
        packed = pack(seq[i:i+len(pattern)], pattern)
        # We count the k-mer as masked if any of its bases are masked.
        masked = packed != packed.upper()
        # Normalize the case of the packed k-mer
        packed = packed.lower()
        if 'n' in packed:
            # Gap/unknown base occurs
            continue
        yield packed, masked

def countSpacedKmers(seq, pattern):
    """Returns a dict from spaced k-mers to values of type:
    [number of occurences, number of masked occurences].

    >>> c = countSpacedKmers('aCTGactg', [True, False, True])
    >>> c == {'at': [2, 2], 'cg': [2, 1], 'ta': [1, 1], 'gc': [1, 1]}
    True
    """
    counts = defaultdict(lambda: [0, 0])
    for kmer, masked in yieldSpacedKmers(seq, pattern):
        counts[kmer][0] += 1
        if masked:
            counts[kmer][1] += 1
    return counts

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
    count = defaultdict(lambda: [0, 0])
    for header, seq in fastaRead(opts.fasta):
        # Add in the counts for this sequence.
        # For speed we should obviously just be passing the count
        # dictionary in rather than adding to it later, but this way is simpler.
        for kmer, v in countSpacedKmers(seq, opts.pattern).iteritems():
            count[kmer][0] += v[0]
            count[kmer][1] += v[1]
    with open('dump', 'w') as f:
        f.write('seed\tnumOccurrences\tnumMaskedOccurrences\n')
        for key, val in count.iteritems():
            occurrences = val[0]
            maskedOccurrences = val[1]
            f.write('%s\t%s\t%s\n' % (key, occurrences, maskedOccurrences))

if __name__ == '__main__':
    main()
