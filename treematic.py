#! /usr/bin/env python

# File: treematic.py Author: Dylan Schwilk (www.pricklysoft.org)

#  $Date: 2008-04-28$

#  Copyright (C) 2008 Dylan W. Schwilk www.pricklysoft.org,
#  Algorithms orginally by Cam Webb

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

"""Python implementation of Cam Webb's phylomatic application.  Variation allows non species level tips.

    citation needed!
    
"""

__version__ =    '''1.1'''
__program__ =    '''treematic.py'''
__author__  =    '''Dylan Schwilk (www.pricklysoft.org)'''
__usage__   =    '''treematic.py [options] [tree_file]'''


import newick
from phylotree import PhyloTree
import logging
logging.basicConfig(format='%(levelname)s: %(message)s')
phylo_logger = logging.getLogger('phylo_logger')


def phylomatic(mtree,taxa, normalize=False, prune=True):
    '''Run phylomatic replacement algorithm on megatree (mtree) using taxa
    list.'''
    newtree = mtree.copy()
    for taxon in taxa:
        matched = False
        for i, lab in enumerate(taxon):  # each taxon list must already be in order species, genus, family, etc
            for node in newtree.postorder():
                if node.label == lab :
                    matched = True
                    phylo_logger.info("FOUND MATCH: " + node.label)
                    toadd = taxon[0:i]
                    toadd.reverse()
                    n = node
                    for t in toadd:
                        c =  PhyloTree(label=t)
                        n.add_child(c)
                        #print "ADDED ", c.label, " TO ", n.label
                        n = c
                    break
            if matched : break
        if not matched : phylo_logger.warning("NO MATCH: " + "/".join(taxon))

    if prune :
        termlist = map(lambda x :x[0], taxa)
        newtree.prune_to_taxa(termlist,normalize=normalize)
    elif normalize :
        newtree.normalize()

    return newtree

def read_taxa_file(f):
    """Read taxa file in phylomatic format. Returns list of lists of form:
           [ [species, genus, family], [species2,genus2, family2] ]"""
    taxa = []
    try :
        for l in open(f).readlines():
            l = l.strip()
            l = l.split("/")
            map(lambda x:x.strip(),l)
            l.reverse() # so species is first
            taxa.append(l)
    except IOError :
        phylo_logger.error('Error reading taxa file, %s' % args[0])
    return taxa


def main():
    '''Command line program.  '''
    import sys   
    from optparse import OptionParser

    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
    parser.add_option("-t", "--taxa", action="store", type="string", \
                      dest="taxa_file",  default = '', help="taxa file in format: family/genus/species")
    parser.add_option("-n", "--normalize", action="store_true", dest="normalize", default=False,
					  help="Normalize resulting tree (collapse nodes with single child), default=%default")
    parser.add_option("-f", "--full", action="store_false", dest="prune", default=True,
					  help="Output full mega tree (do not prune to taxa list)")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")    

    (options, args) = parser.parse_args()

    if options.verbose:
        phylo_logger.setLevel(logging.INFO)
    
    if len(args) == 1 :
        try :
            src = open(args[0]).read()
        except IOError:
            phylo_logger.error('Error reading file, %s' % args[0])
            sys.exit()
    else :
        src = sys.stdin.read()

    trees = newick.read_trees(src)

    # now get taxa
    if options.taxa_file:
        taxa = read_taxa_file(options.taxa_file)

            

    #for n   in trees[0].postorder(): print n.ulabel()

    result = []
    for tree in trees:
        result.append(phylomatic(tree,taxa,options.normalize,options.prune))

    for tree in result:
        #pass
        print "%s;" % tree
        #tree.normalize()
        #print tree.children[0].label
       #     if len(node.children)<2 : print node.label
#        print tree.leaves()
    return 0



if __name__== "__main__":
    main()
