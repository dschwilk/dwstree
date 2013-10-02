#! /usr/bin/env python

# File: branch_lengths.py
# Author: Dylan Schwilk (www.pricklysoft.org) and Helene Morlon
# 2008-09-18

# Copyright (C) 2003, 2004, 2005, 2008 Dylan W. Schwilk and Helene Morlon

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


"""Module for branch length manipulations

    Tree argument to all functions must provide the interface of a PhyloTree.

    TODO: list functions and usage.
    
"""

__version__ = "3.0"
__program__ =    '''branch_lengths.py'''
__author__  =    ['''Dylan Schwilk''','''Helene Morlon''']
__usage__   =    '''branch_lengths.py [options] [tree_file]'''

import numpy

import newick
import logging, random, math
phylo_logger = logging.getLogger('phylo_logger')

### utility functions
def _rand_truncated_exp(lambd, x0):
    """ Random number generator: returns truncated exponential deviate with
    mean at 1.0 /lambd and truncated at x0"""
    return -1.0/lambd*math.log(1.0 - random.random() * (1.0-math.exp(-lambd*x0)))

## Node age distribution functions
## These functions are passed to bladj smoothing function
def node_age_bladj_original(start,stop,n):
    """Returns node ages by original (non-stochastic) bladj algorithm"""
    dist = stop-start
    return [start+(i+1)*dist/(n+1) for i in range(n)]

def node_age_uniform(start,stop, n):
    """Returns n node ages (breakpoints on line start -- stop)."""
    return sorted([random.uniform(start,stop) for i in range(n)])

def node_age_exponential(start, stop, n, alpha, reverse=False):
    """Returns n node ages (breakpoints on line start -- stop) according to a
    truncated exponential distribution with mean = alpha*(stop-start), so lamda
    = (1/alpha)(stop-start). If reverse is true, origin of exponential pdf is
    ancestor age rather than descendant."""
    dist = stop-start
    #lambd = dist*1.0/alpha
    if reverse :
        return sorted([stop - _rand_truncated_exp(1.0/(alpha*dist),dist) for i in range(n)])
    else :
        return sorted([start + _rand_truncated_exp(1.0/(alpha*dist),dist) for i in range(n)])

## public functions
def bl_ages(tree, ages=None):
    """"Transform ages stored in ages dict to branch lengths. Dictionary must
    have a key for every node ulabel (including tips)"""
    if ages:
        for node in tree.postorder():
            if not node.is_tip():
                for child in node.children:
                    child.bl = ages[node.ulabel()] - ages[child.ulabel()]
    else :
        for node in tree.postorder():
            if not node.is_tip():
                for child in node.children:
                    child.bl = node.age - child.age
    
def bl_one(tree):
    """Set all branch lengths to one. Root set to zero"""
    for node in tree :
        if node.is_root():
            node.bl=0.0
        else :
            node.bl=1.0

def bl_grafen(tree):
    """Set branch lengths according to Grafen method where
    node height is proportional to number of descendants. Root set to zero.
   """
   # total_height = len(tree.descendants())  # Should I normalize?
    for node in tree:
        if node.is_root() :
            node.bl=0
        else :
            node.bl = float(len(node.parent.descendants()) - len(node.descendants())) #/ total_height

def bl_topo(tree):
    """Sets branch lengths according to topo method (macclade view)
    where node height is proportional to the maximum number of nodes
    descendant. Root bl set to zero."""
    height = lambda n : max(n.nodes_to_tips())
  #  total_height = height(tree)
    for node in tree :
        if node.is_root() :
            node.bl=0
        else :
            node.bl = float(height(node.parent) - height(node)) #/ total_height

# def bl_bladj(tree, the_age_dict, age_dist_func = node_age_bladj_original, all_tips_zero=False, assume_node_labels=False):
#     """Set node ages according to bladj algorithm and fixed nodes. Root age
#     must be fixed in ages_dict. This version allows using bladj on trees with
#     tips that are not species, i.e allows for assigning branch lengths before
#     pruning trees. This function will modify the tree topology if all_tips_zero
#     is False, not merely branch lengths because it will delete entire branches
#     that lack any fixed node age information. Based on the description of Cam
#     Webb's bladj algorithm in phylocom and email discussions with Cam.

#     Notes on version2: Function now accepts an "age distribution function" to
#     set unfixed node ages between two known age nodes. The default,
#     `_node_age_bladj_original` behaves as original: ages are evenly spaced
#     between known ages. Other functions can produce other distributions (such as
#     stochastic uniform or truncated exponential)

#     Function modifies tree, but not age_dict
#     """
#     age_dict = the_age_dict.copy()  # needed to avoid problem with changing orginal? How to eliminate this copy?

#     if age_dict.get(tree.ulabel(), 0) == 0 :  raise(ValueError("Error in bladj:  root node must have a fixed age"))

#     if all_tips_zero:  ## assume all tips are age=0, regardless of naming convention
#         for t in tree.leaves():
#             age_dict[t.ulabel()]=0
#     else : # set only TRUE species terminal taxa to age=0. Terminal taxa recognized by
#            # having an underscore in the node label.
#         for t in tree.leaves():
#             if t.label:  ## don't use ulabel()!
#                 if "_" in t.label :
#                     age_dict[t.ulabel()]=0

#     _prune_unaged_tips(tree, age_dict)

#     # main loop
#     for node in tree :
#        if not (node.ulabel() in age_dict) : # Only work on non fixed nodes
#            assert(not node.is_tip())  # Should never have unaged tip.
#            parent_age = age_dict[node.parent.ulabel()]
#            fd = _fixed_descendants(node,age_dict) # list of all line-of-sight descendents with fixed age
#            mindist = parent_age - age_dict[fd[0].ulabel()]
#            minindex = 0

#            # Find descendent with fixed age nearest in age to this node's parent
#            for j,d in enumerate(fd[1:]):
#               i=j+1
#               dist = parent_age - age_dict[d.ulabel()]
#               if dist == mindist :  # We've found equivalent fixed descendents
#                   if len(_nodes_to_fixed_parent(d.parent, age_dict)) \
#                           < len(_nodes_to_fixed_parent(fd[minindex].parent, age_dict)):
#                       minindex=i    # we choose shortest intervening node distance if age dist is same
#               elif dist < mindist:  
#                   mindist = dist
#                   minindex = i
#            desc_age = age_dict[fd[minindex].ulabel()]  # winning node, closest fixed age
#            tofix = _nodes_to_fixed_parent(fd[minindex].parent, age_dict)

#            # set node ages:
#            new_ages = age_dist_func(desc_age,parent_age,len(tofix))
#            for i,f in enumerate(tofix):
#                age_dict[f.ulabel()] = new_ages[i]

#     bl_ages(tree, age_dict)  # assign the new branch lengths via the node age method

# def _add_ages_from_dict(tree,age_dict):
#     for node in tree:
#         node.age = age_dict.get(node.label,None)
#         print node.label, node.age
# #        if node.label : node.age = age_dict[node.label]
# #        else : node.age = None

def bl_bladj(tree, age_dict, age_dist_func = node_age_bladj_original, all_tips_zero=False):
    """Set node ages according to bladj algorithm and fixed nodes. Root age
    must be fixed in ages_dict. This version allows using bladj on trees with
    tips that are not species, i.e allows for assigning branch lengths before
    pruning trees. This function will modify the tree topology if all_tips_zero
    is False, not merely branch lengths because it will delete entire branches
    that lack any fixed node age information. Based on the description of Cam
    Webb's bladj algorithm in phylocom and email discussions with Cam.

    Notes on version2: Function now accepts an "age distribution function" to
    set unfixed node ages between two known age nodes. The default,
    `_node_age_bladj_original` behaves as original: ages are evenly spaced
    between known ages. Other functions can produce other distributions (such as
    stochastic uniform or truncated exponential)

    Function modifies tree, but not age_dict. Does not return a value.
    """
    # assign known ages .. For now, store temp ages in .age data member of each node.
    for node in tree:
        if age_dict.has_key(node.label):
            node.age = age_dict[node.label]
        elif node.is_tip() :
            if all_tips_zero or "_" in node.label:
                node.age = 0.0
            else :
                node.age = None
        else :
            node.age = None

    if tree.age is None  :  raise(ValueError("Error in bladj:  root node must have a fixed age"))
    _prune_unaged_tips(tree)  ## get rid of unageable regions of tree

    # main loop to find unaged nodes
    for node in tree :
       if node.age is None : # only work on non fixed nodes
           #assert(not node.is_tip())  # Should never have unaged tip.
           
           # Find descendent with fixed age nearest in age to this node's parent           
           fd = _fixed_descendants(node) # list of all line-of-sight descendents with fixed age
           
           mindist = node.parent.age - fd[0].age 
           minindex = 0
           i = 1
           while i < len(fd):
              dist = node.parent.age - fd[i].age
              if dist == mindist :  # We've found equivalent fixed descendents
                  if len(_nodes_to_fixed_parent(fd[i].parent)) \
                          < len(_nodes_to_fixed_parent(fd[minindex])):
                      minindex=i    # we choose shortest intervening node distance if age dist is same
              elif dist < mindist:  
                  mindist = dist
                  minindex = i
              i += 1
              
           desc_age = fd[minindex].age  # winning node, closest fixed age
           tofix = _nodes_to_fixed_parent(fd[minindex])

           # set node ages:
           new_ages = age_dist_func(desc_age,node.parent.age,len(tofix))
           for i,f in enumerate(tofix):
               f.age = new_ages[i]

    bl_ages(tree)  # assign the new branch lengths via the node age method
    
###############################################################
# private functions

def get_age_dict(filename):
    """Read fixed ages from text"""
    f = open(filename)
    a_dict= {}
    ages = f.readlines()
    f.close()
    for line in ages :
        clade, age = line.split()
        clade = clade.strip()
        age = float(age.strip())
        a_dict[clade] = age
    return a_dict

# def _prune_unaged_tips(tree, age_dict, stop=False):
#     """Remove all lineages with no fixed ages"""
#     prunelist = map(lambda y: y.label , filter(lambda x: not age_dict.has_key(x.label), tree.leaves()))
#     #print prunelist
#     tree.prune_taxa(prunelist)
    
def _prune_unaged_tips(tree):
    """Remove all lineages with no fixed ages"""
    prunelist = map(lambda y: y.label , filter(lambda x: x.age is None, tree.leaves()))
    #print prunelist
    tree.prune_taxa(prunelist)

    
# def _fixed_descendants(node, age_dict, vect=None):
#     """ Returns list of descendants that have ages fixed in age_dict """
#     if vect == None: vect = []
#     if (node.ulabel() in age_dict):
#         vect.append(node)
#         return vect
#     else:
#         for child in node.children:
#             _fixed_descendants(child,age_dict,vect)
#     return vect

def _fixed_descendants(node, vect=None):
    """ Returns list of descendants that have ages fixed in age_dict """
    if vect == None: vect = []
    if not node.age is None:
        vect.append(node)
        return vect
    else:
        for child in node.children:
            _fixed_descendants(child,vect)
    return vect

# def _nodes_to_fixed_parent(node, age_dict, vect=None):
#     """return a vector containing all ancestor nodes excluding most recent ancestor with fixed age
#     """
#     if vect == None: vect = []
#     if age_dict.has_key(node.ulabel()):
#         return vect
#     else :
#         vect.append(node)
#         _nodes_to_fixed_parent(node.parent, age_dict, vect)
#     return vect


def _nodes_to_fixed_parent(node):
    """return a vector containing all ancestor nodes excluding most recent ancestor with fixed age
    """
    vect = []
    par = node.parent
    while  par.age is None :
        vect.append(par)
        par = par.parent
    return vect

############################################################################
# Command line script
############################################################################

def main():
    '''Command line program. Main argument to script is a tree file in newick
    format (or in NEXUS format if the -n option is used). The script performs one of the
    branch length algorithms on all the trees in that file and prints the changed
    file.'''
    import sys   
    from optparse import OptionParser
    logging.basicConfig()

    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
    parser.add_option("-n", "--nexus", action="store_true", dest="nexus", \
                      default = 0 , help="Read trees from NEXUS file rather than newick tree file")  
    parser.add_option("-t", "--topo", action="store_true", dest="topo", \
                      default = 0 , help="Set branch lengths by MacClade method")    
    parser.add_option("-o", "--one", action="store_true", dest="one", \
                      default = 0 , help="Set branch lengths to one")  
    parser.add_option("-g", "--grafen", action="store_true", dest="grafen", \
                      default = 0 , help="Set branch lengths by grafen method") 
    parser.add_option("-s", "--smooth", action="store", type="string", \
                      dest="ages_file",  default = '', \
                      help="Smooth ages between fixed nodes according to bladj algorithm using ages in file")
    parser.add_option("-u", "--uniform", action="store_true", dest="uniform_age_dist",  default = 0, \
                      help="Distribute unknown ages by random uniform distribution. requires -s option.")
    parser.add_option("-e", "--exponential", action="store", dest="exp_age_dist",  type="float", default = 0, \
                      help="Distribute unknown ages by truncated exponential distribution with a mean of (EXP_AGE_DIST*distance).  Negative indicates that skewness should be reversed.  Requires -s option")
    parser.add_option("-r", "--replicates", action="store", dest="nreps",  type="int", default = 1, \
                      help="Number of replicate trees to output (only sensible with -u or -e options)")
    parser.add_option("-a", "--ages", action="store_true", dest="ageout",  \
                      default = 0, help="Output node ages")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")    

    (options, args) = parser.parse_args()

    if options.verbose:
        phylo_logger.setLevel(logging.INFO)
    
    if len(args) == 1 :
        try :
            src = open(args[0]).read()
        except IOError:
            phylo_logger.error("Error reading tree file, %s" % args[0])
            sys.exit()
    else :
        src = sys.stdin.read()

    if options.nexus:
        from nexus_doc import NexusDoc  # no reason 
        nxdoc = NexusDoc(log = None)
        nxdoc.load(src)
        trees = nxdoc.Trees()
    else :
        trees = newick.read_trees(src)

    result_trees= []
    
    for tree in trees:
        if options.one :
            phylo_logger.info("Setting branch lengths to one")
            bl_one(tree)
        elif options.grafen :
            phylo_logger.info("Setting branch lengths by Grafen method")
            bl_grafen(tree)
        elif options.topo :           
            bl_topo(tree)
            phylo_logger.info("Setting branch lengths by minimal extension method")
        elif options.ages_file :
            try:
                a_dict = get_age_dict(options.ages_file)
            except IOError:
                phylo_logger.error( "Error reading ages file: %s" % options.ages_file)
                sys.exit()
            phylo_logger.info("Running bl_bladj")
            age_dist_func = _node_age_bladj_original  # default
            if options.uniform_age_dist:
                phylo_logger.info("Running bl_bladj with uniform age distribution")
                age_dist_func = _node_age_uniform
            if options.exp_age_dist != 0:
                if options.exp_age_dist < 0 :
                    reverse = True
                    options.exp_age_dist = abs(options.exp_age_dist)
                else :
                    reverse = False
                phylo_logger.info("Running bl_bladj with exponential age distribution, alpha=%f" % options.exp_age_dist)
                def node_age(start,stop, n):
                    return _node_age_exponential(start,stop,n,options.exp_age_dist, reverse)
                age_dist_func = node_age
            for i in range(options.nreps):
                newtree = tree.copy()
                bl_bladj(newtree, a_dict, age_dist_func, False)
                result_trees.append(newtree)
        elif options.ageout :
            phylo_logger.info("Printing node ages")
            ages = tree.node_ages()
            for (l,a) in ages:
                print l,a

    if options.ageout:
        return 0
    
    if options.nexus:
        trees = trees + result_trees
        print nxdoc  
    elif len(result_trees) > 0:
        for t in result_trees:
            print t, ";"
    else:
        for tree in trees:
            print tree, ";"
    return 0

        
def test():
    #TODO: add tests for bladj
    
    "Test functions"
    from newick import create_tree
    
    treestr = "((((incanus:5.468750,cordulatus:5.468750):5.468750,(((leucodermis:3.125000,cyaneus:3.125000):3.125000,thyrsiflorus:6.250000):3.125000,((arboreus:3.906250,spinosus:3.906250):3.906250,(integerrimus:6.250000,(hearstiorum:4.687500,(tomentosus:3.125000,(oliganthus:1.562500,impressus:1.562500):1.562500):1.562500):1.562500):1.562500):1.562500):1.562500):1.562500,((papillosus:4.687500,foliosus:4.687500):4.687500,(diversifolius:6.250000,(griseus:3.125000,lemmonii:3.125000):3.125000):3.125000):3.125000)euceanothus:16.500000,(((((((divergens:0.642857,purpureus:0.642857):0.642857,prostratus:1.285714):0.642857,megacarpus:1.928571):0.642857,maritimus:2.571429):0.642857,masonii:3.214286):0.642857,cuneatus:3.857143):0.642857,crassifolius:4.500000)cerastes:24.500000)ceanothus:1.000000"

    tree = create_tree(treestr)
    tree2 = tree.copy()
    bl_topo(tree)
    ages = tree.node_ages()
    a = {}
    for i, node in enumerate(tree2.postorder()) :
        a[node.ulabel()] = ages[i]
    bl_ages(tree2, a)
    print tree.write(True)
    print tree2.write(True)
    if tree.write(True)==tree2.write(True) :
        print "Success"

    ltree =  create_tree(treestr)
    node_ages = {'ceanothus': 29.0, 'cerastes': 4.5, 'euceanothus' : 12.5}
    bl_bladj(ltree, node_ages)
    
    print ltree

if __name__== "__main__":
    main()
    
