#! /usr/bin/env python

# File: explore-bl.py
# Author: Dylan Schwilk
# Copyright (C) 2009 Dylan W. Schwilk

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
rarefaction.py. Script to do phylogenetic distance by taxonomic diversity
rarefaction curves. Command line options exist for running. Also check
constants defined below (eg root of tree).
"""

__version__ = "1"
__needs__ = '2.4'
__author__ = "Dylan W. Schwilk"
__program__ =    '''rarefaction.py'''
__usage__   =    '''rarefaction.py [options] [tree_file]'''

import phylotree, newick
import random,logging, math
import numpy
from parallel import Parallel, delayed
#from scipy import linspace, polyval, polyfit, sqrt, stats, randn
phylo_logger = logging.getLogger('phylo_logger')

### parameters
MAX_RARE_POINTS=32
LOG_S_VALUES=10.0**numpy.linspace(0.25,4,MAX_RARE_POINTS)  # to get 32 points equally spaced in log space
S_VALS=[int(i) for i in LOG_S_VALUES]

# print logspace
NREPS=100  ## number of random td-pd samples per richness

##params for bl exponential explore
ALPHAS=[0.2,0.4,0.6,0.8]
#ALPHAS=[0.2,0.8]

## defaults  -- may be changed by command line scripts
##params for random resolution explore
TOPO_REPS=1000# n of randomized topologies (resolutions, bls, etc)
REROOT_TREE="angiosperms" ## common root for all analyses
FIXED_INTERCEPT=179.0  # temporary

############################################################################
# Command line script
############################################################################

def linear_reg_old(x,y):
    """Linear regression. x and y are numpy arrays"""
    (ar,br)=numpy.polyfit(x,y,1)
    return(ar,br)

def linear_reg(x, y):
    """ returns slope, intercept and r-squared as tuple"""
    coeffs = numpy.polyfit(x, y, 1)

    # r-squared
    p = numpy.poly1d(coeffs)
    # fit values, and mean
    yhat = [p(z) for z in x]
    ybar = sum(y)/len(y)
    ssreg = sum([ (yihat - ybar)**2 for yihat in yhat])
    sstot = sum([ (yi - ybar)**2 for yi in y])
    rsquared = ssreg / sstot
    return (coeffs[0],coeffs[1],rsquared)

def linear_reg_fixed_intercept(x,y,y0):
    """Calculate least squares slope for regression with fixed y intercept =
    y0. x and y should be numpy arrays, y0 should be a scalar. Simple version
    below probably faster."""
    xm = numpy.vstack([x, numpy.zeros(len(x))]).T  # set up coefficient matrix for numpy.linalg.lstsq
    ynew = y-y0
    coeffs,residuals,rank,s = numpy.linalg.lstsq(xm,ynew)
    return coeffs[0]  # just return slope since that is all we estimated

def linear_reg_fixed2(x,y,y0=0.0):
    """Simple slope estimate for fixed regression (regression through origin"""
    ynew = y-y0
    return sum(x*ynew)/sum(x*x)


def tdpd(tree,taxa,nreps):
    """For each value of s, 1 - ntaxa, pull s taxa randomly and calculate pd.
    Repeat nreps times and return numpy vectors for s and pd. There are two
    values of s for which pd is fixed. s=1 (age of tree) and s=len(taxa) (sum
    of branch lengths for tree."""
    points = [i for i in S_VALS if i < len(taxa)]  # use s values spaced evenly in log space
    #points = range(1,len(taxa))  # use all points

    ntaxa = len(taxa)
    if points[-1] != ntaxa: points.append(ntaxa)  ## add full richness value if not in logspaced list
    
    alen = nreps*(len(points) -2 )  # start with array shorter by nreps*2 to allow for s=1 and s=ntax at end
    td = numpy.empty(alen)
    pd = numpy.empty(alen)

    # THis loop below is where these simulations spend the bulk of their time
    i = 0                                     
    for s in points[1:-1]:  #leave off s=1 and s = ntaxa since no need to call random on those
        for r in range(nreps):
            thetax = random.sample(taxa,s)
            td[i] = s
            pd[i] = tree.pd(thetax)
            #print td[i], pd[i]
            i = i+1
    s0vals = tree.pd(random.sample(taxa,1)) * numpy.ones(nreps) ## create nreps values of pd at s=1, assume altrametric
    s1vals = tree.pd(taxa) * numpy.ones(nreps)              ## nreps at s=ntaxa

    td = numpy.concatenate( (td,  numpy.ones(nreps), numpy.ones(nreps)*ntaxa) )
    pd = numpy.concatenate( (pd, s0vals, s1vals) )
    return (td, pd)  

def reroot(tree, lab):
    for node in tree:
        if node.label == lab:
            tree = node
            tree.parent=None
            return tree
    return None # root not found

def read_species_from_taxa_file(fname):
    """Read tips names from a phylomatic/treematic taxa file. Throws away
    genus,family names per taxon."""
    taxa = set()
    try:
        src = open(fname).readlines()
        for l in src :
            l = l.strip()
            l = l.split("/")
            map(lambda x:x.strip(),l)
            taxa.add(l[-1].strip())
    except IOError:
        phylo_logger.error("Error reading taxafile, %s" % options.taxa_file)
        sys.exit()
    return taxa

def getz(t,tx,reps):
    """Just a function to wrap tdpd and following regression together so as to
    have a clean single function to send to parallel processes, etc. Choose
    which regression function to call be commenting out others."""
    x,y = tdpd(t,tx,reps)
#    return linear_reg(numpy.log10(x), numpy.log10(y))
    #return linear_reg_fixed_intercept(numpy.log10(x),numpy.log10(y),numpy.log10(FIXED_INTERCEPT))
    return linear_reg_fixed2(numpy.log10(x),numpy.log10(y),numpy.log10(FIXED_INTERCEPT))


def main():
    """Command line program."""
    import sys   
    from optparse import OptionParser
    from branch_lengths import bl_bladj, get_age_dict, node_age_uniform, node_age_exponential, node_age_bladj_original

    logging.basicConfig()

    parser = OptionParser(usage=__usage__, version ="%prog " + __version__)
    parser.add_option("-a", "--ages-file", action="store", type="string", \
                      dest="ages_file",  default = '', \
                      help="node ages file")
    parser.add_option("-t", "--taxa", action="store", type="string", \
                      dest="taxa_file",  default = '', help="taxa file in format: family/genus/species")
    parser.add_option("-s","--simulation",action="store", type="string", \
                         help="simulation type: NONE | RESOLVE | BLADJ_UNIFORM | BLADJ_EXP.  Default = NONE", dest="sim_type", default="NONE")
    parser.add_option("-j", "--jobs", action="store", type="int", \
                      dest="njobs",  default = 1, \
                      help="Number of jobs for parallel processing. Default = 1")
    parser.add_option("-r", "--reps", action="store", type="int", \
                      dest="toporeps",  default = TOPO_REPS, \
                      help="Number of replicate simulated phylogenies (resolutions or branch lengths) to explore. Default=%d" % TOPO_REPS)        
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
       
    tree = newick.read_trees(src)[0]  ## should contain 1 megatree with full
                                      ## taxa as methods apply branch lengths
                                      ## and resoltuion before pruning to
                                      ## species set in taxa file.
    # read taxa from phylomatic-style taxa file
    taxa = read_species_from_taxa_file(options.taxa_file)
    phylo_logger.info("Read taxa file with %d taxa" % len(taxa))
    
    # read the ages file for bladj smoothing
    age_dict = get_age_dict(options.ages_file)
 
    ## reroot to REROOT_TREE (default=angiosperms)
    tree = reroot(tree, REROOT_TREE)

    ### do sensitivity analysis according to type
    testtrees = []  #vector of alternative phylogenies
    params = []     # vector of parameters for each phylogeny
    if options.sim_type == "NONE":
        bl_bladj(tree, age_dict) # standard bladj
        tree.prune_to_taxa(taxa)
        tree.normalize()
        #print tree
        for i in range(options.toporeps):        
            testtrees.append(tree.copy())
            params.append((options.taxa_file, options.sim_type))
    elif options.sim_type == "RESOLVE":
        for i in range(options.toporeps):
            ntree = tree.copy()
            ntree.resolve() # random resolution
            bl_bladj(ntree, age_dict) # original bladj
            ntree.prune_to_taxa(taxa)
            ntree.normalize() ## reduce nodes to make traversing faster
            testtrees.append(ntree)
            params.append((options.taxa_file, options.sim_type))
    elif options.sim_type == "BLADJ_UNIFORM" :
        for i in range(options.toporeps):
            bl_bladj(tree, age_dict, node_age_uniform)
            tree.prune_to_taxa(taxa)
            tree.normalize() ## reduce nodes to traverse
            testtrees.append(tree.copy())
            params.append((options.taxa_file, options.sim_type))            
    elif options.sim_type == "BLADJ_EXP":
        for alpha in ALPHAS:
            for reverse in [False,True]:
                agefunc = lambda start,stop,n : node_age_exponential(start,stop,n,alpha,reverse)
                for i in range(options.toporeps):
                    bl_bladj(tree, age_dict, agefunc)
                    tree.prune_to_taxa(taxa)
                    tree.normalize() ## reduce nodes to traverse
                    testtrees.append(tree.copy())
                    params.append((options.taxa_file, options.sim_type, alpha, reverse))
    else:
        phylo_logger.error("-s options not recognized.  Possible simulations \
                              types are NONE (default), RESOLVE, BLADJ_UNIFORM, or BLADJ_EXP")

    ### TEST CODE TO OUTPUT RAW TD-PD bivariate data ###
    if True :
        for i, t in enumerate(testtrees):
            x,y = tdpd(t,taxa,NREPS)
            for q in range(len(x)):
                print "TREE%s" % i,
                for p in params[i]: print p,
                print x[q], y[q]
            
        #     t.make_pectinate()
        #     print t.write(True) + ";"
        sys.exit()

    #Now that inputs are created, run tppd using using mulitprocessing module if njobs > 1
    phylo_logger.info("Created %d trees.  Starting rarefaction" % len(testtrees))
    if options.njobs == 1:
        for i,t in enumerate(testtrees):
            for p in params[i]: print p,
            #x,y = tdpd(t,taxa,NREPS)
            #z,b,r = getz(t,taxa,NREPS)
            z = getz(t,taxa,NREPS)
            print z #,b,r            
            sys.stdout.flush()
    else:
        results = Parallel(n_jobs=options.njobs)(delayed(getz)(i,taxa,NREPS) for i in testtrees)
        for  i, r in enumerate(results) :
            for p in params[i] : print p,
            print r #[0], r[1], r[2]
              
    
if __name__ == """__main__""":
    main()

