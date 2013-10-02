# File: parsimony.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2005/03/09 00:59:14 $
# $Revision: 1.2 $
# $Source: /home/schwilk/code/python/cactus-pie/nexus/RCS/parsimony.py $
# Copyright (C) 2003, 2004, 2005 Dylan W. Schwilk (www.pricklysoft.org)

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

"""Module for parsimony state reconstructions

   Functions:
   

"""

__author__  =    '''Dylan Schwilk (www.pricklysoft.org)'''
__version__ =    '''$Revision: 1.2 $'''

from phylotree import PhyloTree
from tree_math import *

def linear_parsimony(tree, char,  matrix):
    """Reconstruct ancester states using linear parsimony

    Not yet finished.

    """

    # downpass
    for node in tree.postorder() :
        if node.isTip() :
            node.data[char] = Range(matrix[(node.label,char)], matrix[node.label])
        else :
            node.data[char] = get_range(matrix[node.label], matrix[(node.label,char)])     
    # now for uppass    
    for node in tree :
        if not node.isTip():
            assert False # not implemented


# def squared_change_parsimony(tree, matrix, char):
#     """Squared-change parsimony"""

#     for node in tree.postorder():
#         if node.isTip() :
#             node.data = matrix[(node.label, char)]
#         else :
#             node.data =                


# Main Test function               
if __name__ == '__main__':
    b = Range((52,170.1))
    a = Range((50.1, 151))
    print a, b
    print get_range(a,b)
  
