#! /usr/bin/env python

# File: dorder.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2008/04/18 19:08:16 $
# $Revision: 1.1 $
# $Source: /home/schwilk/code/python/dwstree/dorder.py,v $
# Copyright (C) 2003, 2004, 2005 Dylan W. Schwilk www.pricklysoft.org

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

"""A module for divergence order and synchronized changes tests."""

__usage__ = """usage: dorder.py [options] char1 char2 test [nexus_filename]
    possible values for 'test' :
        sync       --- synchronized changes (SvS) test
        div-age    --- divergence age test
        print-ages --- simply prints list of contrast magnitudes and divergence ages"""

__version__ =    "1.0"
__author__  =    '''Dylan Schwilk'''
__program__ =    '''dorder'''

from icontrasts import  compute_contrasts, prune_missing_vals
from tree_math import average, weighted_average,  sample_rep
import logging
phylo_logger = logging.getLogger('phylo_logger')

# statistical tests
def sync_change_test(tree, matrix, char1, char2, \
                 nrand = 1000, s_contrasts=False, rand_tips=False):
    """Synchronous change test.  This tests the significance of the
    SvS statistic for two characters. Returns a tuple: SvS, expected
    SvS and P-value."""
    
    c1 =  map(abs, compute_contrasts(tree, matrix, char1, s_contrasts, char1))
    c2 =  map(abs, compute_contrasts(tree, matrix, char2, s_contrasts,char1))
    obs_svs = SvS(c1,c2)

    lower_count = 0
    svs_list = []
    for i in range(nrand):
        if rand_tips :
            rtree = get_randomized_tree(tree)
            r1 = map(abs, compute_contrasts(rtree, matrix, char1, s_contrasts,char1))
            r2 = map(abs, compute_contrasts(rtree, matrix, char2, s_contrasts,char1))
        else :
            r1 = sample_rep(c1)
            r2 = sample_rep(c2)
        rand_svs = SvS(r1,r2)
        svs_list.append(rand_svs)
        if rand_svs <= obs_svs : lower_count = lower_count + 1
    P = lower_count
    exp = sum(svs_list) / float(nrand)
    return obs_svs, exp, float(P) / float(nrand)

def div_age_test(tree, matrix, char1, char2, \
                 nrand = 10000, s_contrasts=False, rand_tips=False):
    """compares weighted average of contrasts age for two characters.
    Results returned are mean1, mean2, observed diff, expected diff, p-value"""
    if len(tree.leaves()) < 3 : return 0,0,0,0

    c1 = map(abs, compute_contrasts(tree, matrix, char1, s_contrasts,char1))
    c2 = map(abs, compute_contrasts(tree, matrix, char2, s_contrasts,char1))
    ages = tree.node_ages()
    m1 = weighted_average(ages, c1)
    m2 = weighted_average(ages,c2)

    # now do sig testing
    rlist = []
    count= 0
    obs_diff = m1-m2
    for i in range(nrand):
        if rand_tips :
            rtree = get_randomized_tree(tree)
            r1 = map(abs, compute_contrasts(rtree, matrix, char1, s_contrasts,char1))
            r2 = map(abs, compute_contrasts(rtree, matrix, char2, s_contrasts,char1))
        else :
            r1 = sample_rep(c1)
            r2 = sample_rep(c2)
            
        rm1 = weighted_average(ages, r1)
        rm2 = weighted_average(ages, r2)
        r_diff = rm1-rm2
        rlist.append(r_diff)
        if obs_diff > r_diff : count = count+1
    exp =  average(rlist)
    p = float(count) / float(nrand)
    return m1,m2,obs_diff,exp,p

##########################################################################
# Other functions
##########################################################################

def divergence_ages(tree, matrix, char1, char2, s_contrasts=False):
    c1 = map(abs,compute_contrasts(tree, matrix, char1, s_contrasts,char1))
    c2 = map(abs,compute_contrasts(tree, matrix, char2, s_contrasts,char1))
    a = tree.node_ages()
    #assert len(a) == len(c1)
    return zip(c1,c2,a)
    #for i in range(len(a)):
    #    print "%f\t%f\t%f" % (c1[i],c2[i],a[i])


# functions for above tests

def SvS(c1, c2):
    """return SvS statistic.
    Algorithm according to description by David Ackerly.
    Input should be absolute values of contrasts"""
    mid_c1 = average(c1)  # average or median or midpoint?
    mid_c2 = average(c2)
    Q1 = 0.0;  Q2 = 0.0 # use floating point
    for i in range(len(c1)):
        if c1[i] >= mid_c1 :
            if c2[i] < mid_c2 : Q1 += 1.0
            elif c2[i] >= mid_c2 : Q2 += 1.0
        else :
            if c2[i] >= mid_c2 : Q1 += 1.0	    
    return Q1 / (Q1+Q2)

def get_randomized_tree(tree):
    """Returns new tree, randomizes taxon labels with replacement.  Old tree is unchanged"""
    result = tree.copy()
    lvs = result.leaves()
    labels = map(lambda n:n.label, lvs)
    rlabels = sample_rep(labels)
    for i, l in enumerate(lvs):
        l.label = rlabels[i]
    return result
        
## Command-line program
## --------------------
def main():
    '''Command line program to read trees and character values from a NEXUS
    file and produce test results.'''
    from nexus_doc import NexusDoc
    from optparse import OptionParser
    import sys    


    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
    parser.add_option("-t", "--trees", action="store", type="string", \
                      dest="trees",  default = 'all', help="trees to include (comma separated)")
#    parser.add_option("-o", "--one-tailed", action="store_true", dest="one_tailed", \
#                      default = 1 , help="Do one-tailed significance tests")
    parser.add_option("-r", "--replicates", action="store", type="int", \
                      dest="NRand", default = 1000, help="Number of randomizations")
    parser.add_option("-s", "--standardized-contrasts", action="store_true", dest="s_contrasts", \
                      default = 0 , help="Standardize contrasts by branch lengths")    
    parser.add_option("-i", "--randomize-tips", action="store_true", dest="rand_tips", \
                      default = 0 , help="Do randomization of tips rather than contrasts") 
    parser.add_option("-v", "--verbose", action="store_true", \
                      dest="verbose",  default = 0, help="verbose output")


    # get options
    (options, args) = parser.parse_args()
    if len(args) == 4 :
        try :
            src = open(args[3]).read()
        except IOError:
            phylo_logger.error('Error in file, %s' % args[0])
            sys.exit()
    elif len(args) == 3 :
        src = sys.stdin.read()
    else :
        p = {}
        p['prog'] = __program__
        print __usage__ % p
        sys.exit()
    char1 = args[0]; char2 = args[1];  test = args[2]
        
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
      
    if options.trees == 'all' :
        treesList = nxdoc.TreeNames()
    else :
        treesList = options.trees.split(',')

    if options.verbose : print "\nTest Results - using characters %s, %s\n" % (char1, char2)

    # Prune trees
    prunedTrees = {}
    # if this code section gets moved to a module, this must be
    # TD[name].copy() so as not to modify original tree :
    #for name in treesList :  prunedTrees[name] = nxdoc.TreeByName(name).copy()
    for name in treesList : prunedTrees[name] = nxdoc.TreeByName(name)
    for name in treesList:
        #branch_lengths.topo_bl(prunedTrees[name])
        prune_missing_vals(prunedTrees[name], [char1,char2], CM)

    
    if test == "sync" :
        # do test on pruned trees
        print "Tree\tSvS\tExp_SvS\tP"
        for name  in treesList :
            tree = prunedTrees[name]
            print "%s\t" % name, 
            print "%f\t%f\t%f" % sync_change_test(tree, CM, char1, char2, \
                                 options.NRand, options.s_contrasts)
    elif test == "div-age" :
        print "Tree\tWMean1\tWMean2\tObsDiff\tExpDiff\tp"
        for name  in treesList :
            tree = prunedTrees[name]
            m1,m2,obs_diff, exp,p = div_age_test(tree, CM, char1, char2, \
                                 options.NRand, options.s_contrasts, options.rand_tips)
            if not options.verbose :
                print "%s\t%f\t%f\t%f\t%f\t%f" % (name,m1,m2,obs_diff, exp,p)

    elif test == "print-ages" :
        print "Tree\t%s\t%s\tnode_age" % (char1, char2)
        for name  in treesList :
            tree = prunedTrees[name]
            results = divergence_ages(tree, CM, char1, char2, options.s_contrasts)
            for line in results :
                print name,
                print "\t%f\t%f\t%f" % line
    else :
        print "Unknown test type. Choose among 'sync', 'div-age', and 'print-ages'"
                 
                              
if __name__ == '__main__':
    main()
