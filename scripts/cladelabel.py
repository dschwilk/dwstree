#! /usr/bin/env python

# File: cladelabel.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# Date: 2005/03/09
# Copyright 2003, 2004, 2005 Dylan W. Schwilk

# GNU
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2 
# of the License, or (at your option) any later version.
#   
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details. 
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.


"""Module for labeling clades

   Functions:
   ----------

       - label_clade(tree, taxa, clade_name) finds the clade defined
       by the taxa in the list and assigns the name clade_name to the
       clade when found.

       - mrca(n1,n2) Finds most resent common ancestor of nodes n1
         and n2

   Command line program:
   

"""
__version__ =    "1.2"
__author__  =    '''Dylan Schwilk (www.pricklysoft.org)'''
__usage__ = '''usage: cladelabel.py clade_label_file [tree_file] '''


from dwstree.phylotree import PhyloTree
import logging
cactus_logger = logging.getLogger('cactus_logger')


def ancestor_list(n):
    """Returns a list of n and its ancestors, starting with the
    current node, (n) and ending with the root of the tree
    """
    result = []
    p=n
    while(1):
        result.append(p)
        if p.is_root() : return result
        p = p.parent

def mrca(n1,n2):
    "Finds most recent common ancestor of n1, n2"
    l1 = ancestor_list(n1)
    p=n2
    #print l1
    while(1):
       # print p
        if p in l1 : return p
        if p.is_root() :  break
        p = p.parent
    return None
    

def label_clade(tree, taxa, clade_name):
    """Find the clade containing all taxa in taxa list and label with
    name"""
    leaves = tree.leaves_by_labels(taxa)
    clade = reduce(lambda a,b : mrca(a,b),leaves)
    if clade :
        clade.label = clade_name

def main():
    '''Command line program.  This program reads a nexus file or newick tree file'''

    from nexus_doc import NexusDoc
    import newick
    import sys   
    from optparse import OptionParser

    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
    
    parser.add_option("-n", "--nexus", action="store_true", dest="nexus", 
                      default = 0 , help="Read trees from NEXUS file rather than newick tree file")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")    

    (options, args) = parser.parse_args()
    if options.verbose:
        cactus_logger.setLevel(logging.INFO)
   
    # Get clade labels
    clade_lines = open(args[0]).read().split('\n\n')

    # Get trees
    if len(args) == 2 :
        try :
            src = open(args[1]).read()
        except:
            print 'Error in tree file, %s' % args[0]
            sys.exit()
    else :
        src = sys.stdin.read()

    if options.nexus:
        nxdoc = NexusDoc(log = None)
        nxdoc.load(src)
        trees = nxdoc.Trees()
    else :
        trees = newick.read_trees(src)
                      
    for line in clade_lines:
        name, taxa = line.split(":")
        name = name.strip()
        taxa = taxa.split()

        for tree in trees:
            label_clade(tree, taxa, name)

    if options.nexus:
        print nxdoc
    else :
        for tree in trees:
            print "%s;" % tree.write(bl=True)

    return 0


# Main Test function               
if __name__ == '__main__':
    main()

    

