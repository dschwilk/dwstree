#! /usr/bin/env python

# File: nexus_doc.py
# Author: Dylan Schwilk
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

"""
    Provides NexusDoc class
"""

import nexus_grammar, sys
from nexus_blocks import TreesBlock, CactusBlock, ContinuousBlock, SetsBlock
from nexus_parser import NexusParser
from copy import copy

class NexusDoc(NexusParser):
    """NEXUS file class for CACTUS documents

        This class provides an interface layer between the underlying nexus_parser
        that contains nexus blocks.  This will allow easier refactoring of underlying
        data representation while maintaining a consistent interface.

        The class has several jobs:
            1. provide brief get/set functions for data contained in blocks
              (ex: treeSet() rather than
              self.blocks['SETS'].objects['TREESET'][self.Properties()['TREESET']] )
            2. Provide default data that a cactus program expects so that
               the cactus client program does not need to worry about
               which blocks are present.

       """
    def __init__(self, log=sys.stdout, use_cactus_block = True):
        NexusParser.__init__(self, log, True)
        self.title = ""
        self.addRecognize('TREES', TreesBlock)
        self.addRecognize('CONTINUOUS', ContinuousBlock)
        self.addRecognize('CACTUS', CactusBlock)
        self.addRecognize('SETS', SetsBlock)

        self.treeset = []
        self.pruneset = []
        self.charset = []
        self.missing =  '?'

    def load(self, input):
        try :
            # see if file-like
            input = input.read()
        except:
            pass
        src = nexus_grammar.parser.parse(input, processor = self)

        # determine if cactus block is present.  If not,
        # get data from continuous block
        if not self.blocks.has_key('CACTUS'):
            if self.blocks.has_key('CONTINUOUS'):
                self.blocks['CACTUS'] = CactusBlock('CACTUS', log=self.log)  # empty cactus block 
                self.blocks['CACTUS'].objects['FORMAT'] = copy(self.blocks['CONTINUOUS'].objects['FORMAT'])
                self.blocks['CACTUS'].objects['CHARLABELS'] = copy(self.blocks['CONTINUOUS'].objects['CHARLABELS'])
                self.blocks['CACTUS'].objects['TAXLABELS'] = copy(self.blocks['CONTINUOUS'].objects['TAXLABELS'])
                self.blocks['CACTUS'].objects['DIMENSIONS'] = copy(self.blocks['CONTINUOUS'].objects['DIMENSIONS'])
                self.blocks['CACTUS'].objects['MATRIX'] = copy(self.blocks['CONTINUOUS'].objects['MATRIX'])
            else :
               self.blocks['CACTUS'] = CactusBlock('CACTUS', log=self.log)  # empty cactus block 
               self.missing =  self.blocks['CACTUS'].getObject('FORMAT','MISSING')

        
        # get sets or assign if they are missing
        if not self.blocks.has_key('SETS') :
            self.blocks['SETS'] = SetsBlock('SETS', log=self.log)
        
        try :        
            self.treeset = self.blocks['SETS'].objects['TREESET'][self.Properties()['TREESET']]
        except KeyError :
            self.blocks['SETS'].addAttribute('TREESET',self.Properties()['TREESET'], self.treeset)
        try:
            self.charset = self.blocks['SETS'].objects['CHARSET'][self.Properties()['CHARSET']]
        except KeyError :
            self.blocks['SETS'].addAttribute('CHARSET',self.Properties()['CHARSET'], self.charset)
        try :
            self.pruneset = self.blocks['SETS'].objects['TAXSET'][self.Properties()['PRUNESET']]
        except KeyError :
           self.blocks['SETS'].addAttribute('TAXSET',self.Properties()['PRUNESET'], self.pruneset)
 
         
    def Taxa(self):
        "Return list of all taxa, or empty list if no Matrix"
        try :
            return self.blocks['CACTUS'].objects['TAXLABELS']
        except KeyError:
            return []

    def Tree(self, i):
        '''Return Tree by number, 1-indexed'''
        try :
            return self.blocks['TREES'].objects['TREES'].items()[i-1][1]
        except IndexError, KeyError :
            return None

    def Trees(self):
        return self.blocks['TREES'].objects['TREES'].values()

    def TreeByName(self, name):
        return self.blocks['TREES'].objects['TREES'][name]

    def TreeNames(self):
        try :
            return self.blocks['TREES'].objects['TREES'].keys()
        except KeyError:
            return []

    def CharMatrix(self):
         return self.blocks['CACTUS'].objects['MATRIX']        

    def Char(self, i):
        "Return Character label for number"
        return self.blocks['CACTUS'].objects['CHARLABELS'][i]
        
    def CharNames(self):
        try :
            return self.blocks['CACTUS'].objects['CHARLABELS']
        except KeyError:
            return []

    def Properties(self):
        return self.blocks['CACTUS'].objects['PROPERTIES']

# Operations
# ----------

