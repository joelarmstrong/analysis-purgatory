#!/usr/bin/env python
from argparse import ArgumentParser
from subprocess import check_output
import itertools
import sys

from sonLib.bioio import fastaRead

def kmers(iterable, k):
    iters = itertools.tee(iterable, k)
    for i, iter in enumerate(iters):
        for _ in xrange(i):
            next(iter, None)
    return itertools.izip(*iters)

def parse_line(line):
    fields = line.split()
    seq = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    mut = fields[3]
    if mut[0] != 'S':
        # not a SNP
        old = ''
        new = ''
    else:
        old = mut[2]
        new = mut[3]
    return seq, start, end, old, new

def main():
    parser = ArgumentParser()
    parser.add_argument('mutationBed')
    parser.add_argument('fasta')
    parser.add_argument('--context', type=int, default=1,
                        help='number of bases of context on each side')
    opts = parser.parse_args()
    k = opts.context * 2 + 1

    genome = dict(fastaRead(opts.fasta))
    with open(opts.mutationBed) as f:
        for lines in kmers(f, k):
            early_exit = False
            for line in lines:
                if len(line) == 0 or line[0] == '#':
                    early_exit = True
                    break
            if early_exit:
                continue
            # Center line
            seq, start, end, old, new = parse_line(lines[opts.context])
            if old.upper() == new.upper():
                continue
            # Check for any overlap with the context of previous
            # lines. If there is overlap, we assume this isn't a valid
            # mutation.
            for prev_line in lines[:opts.context]:
                prev_seq, _, prev_end, _, _ = parse_line(prev_line)
                if prev_seq == seq and prev_end >= start - opts.context:
                    early_exit = True
            # Same, but for the next lines.
            for next_line in lines[opts.context+1:]:
                next_seq, next_start, _, _, _ = parse_line(next_line)
                if next_seq == seq and next_start <= end + opts.context:
                    early_exit = True
            if early_exit:
                continue
            if start < opts.context or start >= len(genome[seq]) - opts.context:
                sys.stderr.write("%s %s %s\n" % (len(context), start, len(genome[seq])))
                continue
            try:
                context = genome[seq][start - opts.context:start + opts.context + 1]
            except Exception as e:
                sys.stderr.write("%s\n" % e)
                continue
            assert len(context) == k
            assert context[opts.context].upper() == new.upper()
            context = context[:opts.context] + old + context[opts.context+1:]
            context = context.upper()
            new = new.upper()
            print '%s\t%s' % (context, new)

if __name__ == '__main__':
    main()
