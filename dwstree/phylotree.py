#! /usr/bin/env python

# File: phylotree.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2006/09/14 04:37:50 $
# $Revision: 1.3 $
# $Source: /home/schwilk/code/python/cactus-pie/nexus/RCS/phylotree.py $
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

"""
PhyloTree module.

Provides n-ary tree data structure with functions to facilitate phylogenetic
analyses. Design emphasizes simplicity over generality and assumes one is
working on rooted trees.

Some ideas for this class (but no code) came from Rick Ree's Mavric
Python package (http://bioinformatics.org/mavric/) but I opted not to
use the Felsenstein linked list approach (FNode), althought such an
approach allows non-recursive pre-order traversal.
"""

__version__ = "1.5"
__needs__ = '2.4'
__program__ =    '''treematic.py'''
__author__  =    '''Dylan Schwilk'''
__usage__   =    '''phylotree.py [options] [tree_file]'''


import copy
import logging
import random
phylo_logger = logging.getLogger('phylo_logger')

######################################################################
# Class: PhyloTree
######################################################################
class PhyloTree:
    """Tree class.

       A tree is simply a handle to a node.  Each node contains:
          - bl: branch length
          - label : the label
          - parent: a ref to the node's parent
          - children: list of node's children
    """
    
    def __init__(self, parent=None, bl = 0.0, label = None):
        self.parent = parent
        self.children = []
        self.label = label
        self.bl = bl

    def __iter__(self):
        """The standard preorder traversal iterator."""
        yield self
        for subtree in self.children :
            for child in subtree :
                yield child

    def __repr__(self) :
        return self.write(True)
     
    def delete(self) :
        if self.parent :
            self.parent.children.remove(self)           
        for subtree in self:
            del subtree
        del self
        
    def postorder(self):
        """Postorder traversal of a tree."""
        for subtree in self.children:
            for child in subtree.postorder():
                yield child
        yield self

    def postorder_list(self):
        """returns postorder list of nodes"""  
        r = []
        for n in self.postorder():
            r.append(n)
        return r

    def preorder_list(self):
        """returns preorder list of nodes"""  
        r = []
        for n in self:
            r.append(n)
        return r
    
    def copy(self) :
        '''Deep copy of tree.'''
        if self.is_tip() :
            return copy.copy(self)
        else :
            newnode = copy.copy(self)
            newnode.children = []
            for subtree in self.children :
                newnode.add_child(subtree.copy())
            return newnode

    def ulabel(self):
        '''Returns the node label if it exists or an entire newick tree string
        (without branch lengths) to uniquely identify a node'''
        if self.label: return self.label
        else : return self.write(False)
        
    def labelize(self):
        ''' Adds an automatic label to any node that lacks a label'''
        n=0
        for node in self:
            if not node.label:
                node.label = "node%d" %n
            n = n+1
                
    def nodes_to_tips(self, vect=None, n=0):
        """return a list of how many internodes are between node and its
        leaves """
        if vect == None: vect = []   
        if self.is_tip():
            vect.append(n)
        else:
            for child in self.children:
                child.nodes_to_tips(vect, n+1)
        return vect

    def length_to_tips(self, vect=None, length=0.0):
        """return a list of total lengths between node and its leaves"""
        if vect == None:
            vect = []
            node_length = 0.0
        else:
            node_length = self.bl #or 1.0

        if self.is_tip():
            vect.append(length+node_length)
        else:
            for child in self.children:
                child.length_to_tips(vect, length+node_length)
        return vect

    def node_ages(self):
        """returns list of node ages for nodes in postorder order. Uses maximum 
        age if tips are non-contemporaneous
        """
        ages = []
        for node in self.postorder():
            L = node.length_to_tips()
            ages.append((node.label,max(L)))
        return ages

    def distance_to_root(self):
        """Returns distance from node to the root"""
        if not self.parent is None:
            return self.parent.distance_to_root() + self.bl
        else:
            return 0
        

    def sum_bl(self, include_root=False) :
        """Sum branch lengths in tree"""
        result = 0
        for node in self:
            result += node.bl
        if include_root :
            return result
        else :
            return result - self.bl

    # def pd2(self, l, include_root=False):
    #     """calculate phylogenetic distance (total bl) for a set of taxa in l"""
    #     inodes = set()
    #     result = 0.0
    #     for node in self.postorder():
    #         if  (node in inodes) or (node.label in l) and not node.parent is None:
    #             result = result + node.bl
    #             inodes.add(node.parent)
    #     if include_root: result = result + self.bl
    #     return result
                
    def pd(self,l):
        """calculate phylogenetic distance (total bl) for a set of taxa in l"""
        inodes = set()
        result = 0.0
        for ch in self.children:
            for node in ch.postorder():
                if  (node in inodes) or (node.label in l):
                    result = result + node.bl
                    inodes.add(node.parent)
        return result
           
    def is_root(self) :
        return self.parent is None

    def is_tip(self) :
        return not self.children

    def is_polytomy(self):
        return len(self.children) > 2

    def add_child(self, child):
        """Add child to parent"""
        child.parent = self
        self.children.append(child)

    def unlink_child(self, child):
        """Unlink a child from parent and delete child."""
        self.children.remove(child)
        del(child)

    def descendants(self):
        """Returns a list of all descendants of this node."""
        d = []
        if not self.is_tip():
            d.append(self)
            for child in self.children:
                if child.is_tip(): d.append(child)
                else: d=d+child.descendants()
        return d

    def leaves(self):
        """Returns a list of leaf nodes that are descendant from this
        node.  Returns a list, is not an iterator to allow modifying tree.
        """
        lvs = []
        if not self.is_tip():
            for child in self.children:
                if child.is_tip(): lvs.append(child)
                else: lvs = lvs+child.leaves()
        else : lvs.append(self)
        return lvs

    def leaves_by_labels(self, labels):
        """returns list of tips that match labels in labels."""
        lvs = []
        for node in self:
            if node.is_tip() and node.label in labels:
                lvs.append(node)
        return lvs       

    def resolve(self):
        """Resolves polytomies to arbitrary zero-length bifurcating branches"""
        for node in self:
            random.shuffle(node.children)
            while len(node.children) > 2 :
                newnode = PhyloTree(bl = 0)
                newnode.add_child(node.children.pop())
                newnode.add_child(node.children.pop())
                node.add_child(newnode)
                
    def di2multi(self,tol=1e-08):
        '''Collapse multichotomies. tol is a numeric value giving the tolerance
        to consider a branch length significantly greater than zero.'''
        for node in self:
            for child in self.children:
                if child.bl < tol:
                    child.collapse()

        
    # def collapse(self) :
    #     "Delete node and add node's children to parent. Does not relabel"
    #     assert not self.is_tip()
    #     assert not self.is_root()
    #     for subtree in self.children :
    #         self.parent.add_child(subtree)
    #     self.parent.bl += self.bl
    #     self.parent.unlink_child(self)
    #     del self

    def normalize(self):
        """Removes all nodes with only one child.  Current version can move root."""
        for node in self.postorder_list():   # can't use generator
            if len(node.children) == 1 :                 
                thechild = node.children[0]
                if node.is_root():  # we can't delete root because that is handle to the whole tree
                    # We need to keep self reference to root so copy child here
                    for gc in thechild.children:
                        node.add_child(gc)# give single child's children to root
                        gc.bl += thechild.bl
                    node.unlink_child(thechild) # delete original single child 
                    
                else :  # add the thechild to its grandparent, delete current node
                    thechild.bl += node.bl  # extend to account for deleting this node
                    node.parent.add_child(thechild)
                    node.parent.unlink_child(node)
        return
    
                
    def prune_taxa(self, l, normalize=False) :
        """Prunes taxa in set or list l from tree. Use a set."""
        for n in self.leaves():
            if n.label in l :
                p = n.parent
                p.unlink_child(n)
                while p.is_tip():
                    p.parent.unlink_child(p)
                    p = p.parent
        if normalize:
            self.normalize()

    def prune_to_taxa(self, l, normalize=False) :
        """Prunes tree leaving taxa in set or list l from tree"""
        for n in self.leaves():
            if n.label not in l :
                p = n.parent
                p.unlink_child(n)
                while p.is_tip():
                    p.parent.unlink_child(p)
                    p = p.parent
        if normalize:
            self.normalize()          


        
    def write(self, bl=False):
        result = ''
        if self.children :
            l = len(self.children)
            result =  "("
            for i in range(l) : 
                if i < l-1 : result = result + self.children[i].write(bl) + ','
                else : result = result + self.children[i].write(bl)
            result = result + ')'
        if self.label :
            result += self.label
        if bl : # and self.parent : # root is allowed to have bl
            result += ':%g' % self.bl
        return result

 
    def write2(self,bl=False):
        if self.children :
            r = ["(", ",".join(map(lambda x: x.write(bl), self.children)),")"]
        else :
            r = []
        if self.label :
            r.append(self.label)
        if bl :
            r.append(":%g" % self.bl)
        return "".join(r)
    
    def write3(self,bl=False):
        if self.children :
            r = "("
            for c in self.children[:-1]:
                r += c.write(bl) +","
            r += self.children[-1].write(bl) + ")"
        else :
            r = ""
        if self.label :
            r+=self.label
        if bl :
            r+=":%g" % self.bl
        return r
    
    def make_pectinate(self):
        """
        Order descendant branches according to their size, so largest
        subtrees are first.  For branches with equal number of
        children, order by label.

        TODO: provide reverse sort. Solved: Easiest to make_pectinate then use
        reverse()
        """
        def sort_func(n1, n2):
            if n1.is_tip() and n2.is_tip():
                #return 0  # could have it ignore labels
                return cmp(n1.label, n2.label)
            else:
                n1lvs = n1.leaves()
                n2lvs = n2.leaves()
                l1 = len(n1lvs)
                l2 = len(n2lvs)
                if l1 > l2:
                    return -1
                elif l1 == l2:
                    #return 0
                    lab1 = map(lambda x:x.label, n1lvs); lab1.sort()
                    lab2 = map(lambda x:x.label, n2lvs); lab2.sort()
                    return cmp(lab1, lab2)
                elif l1 < l2:
                    return 1

        # now the actual sort
        for node in self.postorder_list():
            node.children.sort(sort_func)

    def reverse(self) :
        '''Reverse order of all nodes.'''
        for child in self.children:
            if not child.is_tip() : child.reverse()
        self.children.reverse()


    def relabel_taxa(self, thedict):
        '''relabels tips by translating from dictionary.'''
        for l in self.leaves():
            l.label = thedict.get(l.label,l.label) # default is just to keep old label


## End: PhyloTree class
######################################################################



######################################################################
## Utility functions
######################################################################

def equivalent(a, b, with_bl = False):
    """Tests for equivalent trees.
       This is useful for eliminating topologically equivilant trees,
       for example, in implementing a function such as PAUP's
       'condense trees.'
    """
    pa = a.copy()
    pb = b.copy()
    pa.make_pectinate()
    pb.make_pectinate()
    return __clade_equiv(pa,pb, with_bl)


def count_polytomies(tree):
    count = 0
    for n in tree:
        if len(n.children) > 2 : count = count+1
    return count

            

######################################################################
## Private Utility functions
######################################################################

def __clade_equiv(a,b, with_bl = False):
    '''Assumes both clades are sorted by pectinate. Private.'''
    if a.is_tip() and b.is_tip():
        if a.label != b.label : return False
        else : return True
    else :
        if len(a.children) == len(b.children) :
            for i in range(len(a.children)):
                if not __clade_equiv(a.children[i],b.children[i],with_bl):
                    return False
        else :
            return False
    return True


if __name__ == "__main__":
    '''Command line program.  '''
    import sys
    import newick
    from optparse import OptionParser

    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
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


    for tree in trees:
        p = count_polytomies(tree)
        print "N polytomies", p

        t = len(tree.leaves())
        n = len(tree.descendants())+1

        print "Taxa: %d, Nodes: %d" % (t,n)
        #tree.normalize()
        #print tree.write(True) + ";"
        

    
    



    

