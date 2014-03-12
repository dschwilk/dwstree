


import newick
import phylotree
import sys


f = sys.argv[1]

tree = newick.read_trees(open(f).read())[0]


def isorder(s):
    if "_" in s : return False
    if s[-4:] == "ales": return True

def isfamily(s):
    if "_" in s : return False
    if s[-5:] == "aceae": return True

for c in tree:
    if c.label and isorder(c.label): c.children = []

for c in tree:
    if c.label and isfamily(c.label): c.children = []
    

print tree, ";"
