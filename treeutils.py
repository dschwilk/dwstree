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

"""Tree utilities
    
"""

__version__ =    '''1.0'''
__program__ =    '''treeutils.py'''
__author__  =    '''Dylan Schwilk'''
__usage__   =    '''treematic.py [options] [tree_file]'''


import newick
from phylotree import PhyloTree
import logging
logging.basicConfig(format='\n%(levelname)s:\n%(message)s\n')
phylo_logger = logging.getLogger('phylo_logger')


def read_taxa(src):
    """Read taxa file in phylomatic format"""
    taxa = []
    for l in src :
        l = l.strip()
        l = l.split("/")
        map(lambda x:x.strip(),l)
        l.reverse() # so species is first
        taxa.append(l)
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
    parser.add_option("-p", "--prune", action="store", dest="prune", default='',
					  help="Prune to taxa in file")
    parser.add_option("-r", "--reroot", action="store", dest="reroot", default='',
					  help="Reroot tree at node with label")
    parser.add_option("-v", "--verbose", action="store_true", dest="verbose", default=False,
					  help="Print INFO messages to stdout, default=%default")    

    (options, args) = parser.parse_args()

    if options.verbose:
        phylo_logger.setLevel(logging.INFO)
    
    if len(args) == 1 :
        try :
            src = open(args[0]).read()
        except IOError:
            phylo_logger.error('Error reading tree file, %s' % args[0])
            sys.exit()
    else :
        src = sys.stdin.read()

    trees = newick.read_trees(src)

    # now get taxa
    if options.taxa_file:
        taxa = read_taxa(open(options.taxa_file))
 
    result = []
    for tree in trees:
        if options.prune:
            try :
                ptaxa = read_taxa(open(options.prune).readlines())
            except IOError:
                phylo_logger.error('Error reading file, %s' % options.prune)
                sys.exit()
            
            ptaxa = [i[0] for i in ptaxa]
            tree.prune_to_taxa(ptaxa)
            
        if options.reroot:
            for node in tree:
                if node.label and node.label==options.reroot:
                    tree = node
                    tree.bl=0.0
                    tree.parent = None
                    break
 

        if options.normalize:
            tree.normalize()
        result.append(tree)
        
        #result.append(phylomatic(tree,taxa,options.normalize,options.prune))

    ## Print results
    for tree in result:
        print "%s;" % tree
    return 0



if __name__== "__main__":
    main()
