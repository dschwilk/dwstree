#! /usr/bin/env python

# File: icontrasts.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2005/03/09 00:31:29 $
# $Revision: 1.5 $
# $Source: /home/schwilk/code/python/cactus-pie/nexus/RCS/icontrasts.py $
# Copyright (C) 2003, 2004, 2005 Dylan W. Schwilk
# www.pricklysoft.org

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

"""Module for Felsenstein independent contrasts.
"""

__version__ = "1.5"
__author__  = '''Dylan Schwilk'''


import math
from phylotree import PhyloTree
from nexus_doc import NexusDoc
from tree_math import average


def compute_contrasts(tree, matrix, char, adjusted=True, split_char=None):
    """Compute independent contrasts for the phylogenetic tree.  This
    function returns a lists of contrasts for a single character.  It
    assumes that the matrix passed in contains character values
    accessed by tuple: (taxon, char).  Reconstructed character values
    are stored in the node's data dictionary accessed by char.  values
    in dictionary matrix.  Contrasts are returned as a list containing
    normalized contrast for each node in postorder order.  Polytomies
    are handled as suggested by Pagel (1992) (with split_char used to
    split children into two groups"""
    
    result = []
    for node in tree.postorder():
        node.abl = node.bl # add adjusted branch length attribute to each node  [TODO: Better way?]
        if node.is_tip() :
            node.data[char] = matrix[(node.label, char)] # put char data in data attribute
        elif node.is_polytomy():
            result.append(_do_polytomy(node,char, adjusted, split_char))
        else :
            result.append(_do_felsenstein(node,char, adjusted))
    return result


def contrasts_matrix(treeDict, CharMatrix, treeList, charList, adjusted=True):
    '''Produces a results matrix of contrasts for every char in
    charlist (by char label, every tree in treeList (by name).  Returns a dictionary in
    (tree, char) tuples and the values are lists of contrasts.'''
    results = {}
    for name, tree in treeDict.items():
        for char in charList :
            results[(name, char)] = compute_contrasts(tree,CharMatrix, char, adjusted)
    return results
            
def formatMatrix(matrix, treelist, charlist,with_ages=False):
    """Produces string formatted for output"""
    results = ['TREE\t']
    for c in charlist : results.append('%s\t' %c)
    if with_ages : results.append('%s' % 'NodeAge')
    results.append('\n')
    for name,t in treelist.items() :
        nodelist = [x for x in t.postorder() if not x.is_tip()]
        for i, node in enumerate(nodelist):
                results.append("%s\t" % name )
                for c in charlist :
                   results.append( "%f\t" % matrix[(name,c)][i])
                if with_ages :
                    L = node.length_to_tips()
                    results.append("%f\t" % ( sum(L) / float(len(L)) ) )
                results.append( '\n')
    return ''.join(results)  


##############################################################
## Local functions

def mean(values):
    """Return the arithmetic average of the values."""
    return sum(values) / float(len(values))

def _do_felsenstein(node, char, adjusted=True):
    """Get contrasts for this node.  We assume this is a non-terminal
    node with two children."""

    if adjusted :
        v1,v2 = node.children[0].abl , node.children[1].abl
    else :
        v1 = 1.0; v2 = 1.0
    #extend branch length of this node to create correct variance
    node.abl += (v1 * v2) / (v1 + v2)
    #Save reconstructed value for this node
    node.data[char] = (((1.0/v1)*node.children[0].data[char] \
                        + (1.0/v2)*node.children[1].data[char]))/ (1.0/v1 + 1.0/v2);   
    return (node.children[0].data[char] - node.children[1].data[char]) / math.sqrt(v1 + v2)

 
def _do_polytomy(node, char, adjusted=True, split_char=None):
    """Do Felsenstein contrast with Pagel method for polytomies.  If
    split_char is supplied, the polytomy is split according to this
    character, rather than the one used for contrasts."""
    if not split_char : split_char=char
    sort_by_split_char = lambda n1,n2 : cmp(n1.data[split_char], n2.data[split_char])
    get_char = lambda n : n.data[char]
    get_split_char = lambda n : n.data[split_char]
    children = node.children
    children.sort(sort_by_split_char) 
    N = len(children)
    
    middle = N/2
    if N % 2 != 0 :  # if odd number of children, assign median according to mean of split_char
        mean = average(map(get_split_char, children))
        if children[middle].data[split_char] > mean :
            middle = (N/2) + 1
    

    # get Branch lengths of two groups
    # use branch lengths only if adjusted=True
    if adjusted :
        get_abl = lambda n : n.abl
        v1 = sum(map(get_abl, children[:middle] )) / float(N)
        v2 = sum(map(get_abl, children[middle:] )) / float(N)
    else :
        v1 = 1; v2 =1

    node.abl += (v1 * v2) / (v1 + v2)
                   
    # get char values for two groups                   
    c1 = sum(map(get_char, children[:middle])) / float(N)
    c2 =  sum(map(get_char, children[middle:])) / float(N)
         
    node.data[char] = (((1.0/v1)*c1 + (1.0/v2)*c2))/ (1.0/v1 + 1.0/v2);   
    return (c1 - c2) / math.sqrt(v1 + v2)
     
    
def prune_missing_vals(tree, charList, charMatrix)  :
    """Prune all taxa with missing character data"""
    prune_list = []
    taxa = map(lambda l : l.label, tree.leaves())
    for c in charList :
        for taxon in taxa :
            if not charMatrix.has_key((taxon,c)) : prune_list.append(taxon)
    tree.prune_taxa(prune_list)

#def _clean(tree) :
#    """Delete the extra node attribute abl"""
#    for node in tree : del node.abl

## Command-line program
## --------------------
def main():
    '''Command line program to read trees and character values from a
    NEXUS file and produce independent contrasts.'''
    from nexus_doc import NexusDoc
    import math, sys
    
    try:
    # Use 2.3 optparse module
        from optparse import OptionParser
    except (ImportError, AttributeError):
        try:
            from optik import OptionParser
        except (ImportError, AttributeError):
            print """Needs python 2.3 or Greg Ward's optik module."""

    
    usage = "usage: %prog [options] filename"
    parser = OptionParser(usage=usage, version ="%prog " + __version__)
    parser.add_option("-c", "--characters", action="store", type="string", \
                      dest="chars", default = 'all', help="characters to include (comma separated)")
    parser.add_option("-t", "--trees", action="store", type="string", \
                      dest="trees",  default = 'all', help="trees to include (comma separated)")
    parser.add_option("-a", "--age", action="store_true", \
                      dest="age",  default = 0, help="report node ages")
    parser.add_option("-v", "--verbose", action="store_true", \
                      dest="verbose",  default = 0, help="verbose output")



    # get options
    (options, args) = parser.parse_args()
    if len(args) == 1 :
        try :
            src = open(args[0]).read()
        except:
            print 'Error in file, %s' % args[0]
    else :
        src = sys.stdin.read()

    # verbose output
    if options.verbose :
        log = sys.stdout
    else :
        log = None

    nxdoc = NexusDoc(log = log)
    nxdoc.load(src)
    CM = nxdoc.CharMatrix()
    taxa = nxdoc.Taxa()
    characters = nxdoc.CharNames()
      
    if options.chars == 'all' :
        charList = characters
    else :
        charList = options.chars.split(',')

    if options.trees == 'all' :
        treesList = nxdoc.TreeNames()
    else :
        treesList = options.trees.split(',')

    # Prune trees
    prunedTrees = {}
    # if this gets moved to a module, this must be TD[name].copy() so as not to modify original tree :
    # for name in treesList :  prunedTrees[name] = TD[name].copy()
    for name in treesList :  prunedTrees[name] = nxdoc.TreeByName(name)
    for name in treesList: prune_missing_vals(prunedTrees[name], charList, CM)
    
    # do contrasts on pruned trees
    if options.verbose : print "\nContrasts Matrix:\n"
    contrasts = contrasts_matrix(prunedTrees, CM, treesList, charList)
    print formatMatrix(contrasts, prunedTrees, charList, with_ages = options.age)

                          
if __name__ == '__main__':
    main()
