"""Rotate a newick tree to put the leaf with a given label first."""
from argparse import ArgumentParser
import newick

def rotate(node, label):
    if node.name == label:
        return True
    if len(node.descendants) == 0:
        return False
    assert len(node.descendants) == 2
    if rotate(node.descendants[1], label):
        node.descendants = reversed(node.descendants)
        return True
    if rotate(node.descendants[0], label):
        return True

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('newick_file')
    opts = parser.parse_args()
    tree = newick.read(opts.newick_file)[0]
    rotate(tree)
    print(newick.dumps(tree))
