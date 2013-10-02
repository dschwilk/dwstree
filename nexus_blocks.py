#! /usr/bin/env python

# File: nexus_blocks.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2005/03/09 00:44:11 $
# $Revision: 1.2 $
# $Source: C:\\local\\dylan\\code\\python\\cactus-pie\\nexus\\RCS/nexus_blocks.py $
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


__version__ = "0.1"
__needs__ = "2.2"
__author__ = "Dylan W. Schwilk"

"""nexus_blocks module

   Provides several nexus block classes:
    TreesBlock
    ContinuousBlock
    SetsBlock.
"""

from nexus_parser import BlockProcessor, NexusError, InputError, nx_string
import newick, sys, operator
from nexus_dict import NexusDict

##############################################################################
## ContinuousBlock
##############################################################################
class ContinuousBlock(BlockProcessor):
    '''CONTINUOUS block.

       stores the following data:
        charlabels: list of character labels
        matrix: a dictionary of character data
    '''
    def __init__(self, n = 'CONTINUOUS', log = sys.stdout):
        BlockProcessor.__init__(self, n, log)
        self.objects['MATRIX'] = {}  # use regular dictionary here  must look up taxa by case sensitive name
        # defaults
        self.objects['FORMAT'] = NexusDict()
        self.objects['FORMAT']['MISSING'] = '?'
        self.objects['FORMAT']['ITEMS'] = ['AVERAGE']

    def __repr__(self):
        return self.asString()

    def asString(self):    
        from nexus_parser import nx_string
        result = ["BEGIN %s;\n" % self.blockname,]
        result.append(self.writeDimensions())
        result.append(self.writeFormat())
        result.append(self.writeCharlabels())
        result.append(self.writeMatrix())
        result.append('END;\n')
        return ''.join(result)

    def writeDimensions(self):
        result = []
        result.append('DIMENSIONS ')
        for dim, num in self.objects['DIMENSIONS'].items():
            result.append("%s = %d " %(dim, num))
        result.append(";\n")
        return ''.join(result)

    def writeFormat(self):
        result = ["FORMAT ",]
        for name, desc in self.objects['FORMAT'].items():
            result.append("%s = %s " %(name, nx_string(desc)))
        result.append(";\n")
        return ''.join(result)

    def writeCharlabels(self):
        result = ['CHARLABELS',]
        for i, char in enumerate(self.objects['CHARLABELS']):
            result.append(" %s [%d]" % (nx_string(char), i))
        result.append(';\n')
        return ''.join(result)

    def writeMatrix(self):
        result = ['MATRIX\n',]
        for taxon in self.objects['TAXLABELS']:
            result.append("\t%s\t" % nx_string(taxon))
            for char in self.objects['CHARLABELS']:
                value = self.objects['MATRIX'].get((taxon,char), self.objects['FORMAT']['MISSING'] )
                if type(value) == type([]) :
                    result.append('(')
                    for i in value :
                        result.append(str(i))
                    result.append(') ')
                else :
                    result.append("%s\t" % str(value))  # number of missing_val
            result.append('\n')
        result.append(';\n')
        return ''.join(result)

    def __checkCharlabels(self, nc):
        if not self.objects.has_key('CHARLABELS'):
            self.objects['CHARLABELS'] = []
            for i in range(1,nc+1):
                self.objects['CHARLABELS'].append("CHAR%d" % i)
            

    def assignObject(self, category, object_name, object_description, format = None, has_star = 0) :
        '''Handle assignments of type: CATEGORY * NAME = (FORMAT) = description.'''
        raise InputError( category, 'Unrecognized object category')

        
    # commands
    def doCHARLABELS(self, token_list):
        '''Assign character labels'''
        self.objects["CHARLABELS"] = token_list 
        #print token_list

            
    def doMATRIX(self, token_list):
        '''Reads MATRIX command.'''
        
        nc = self.objects['DIMENSIONS']['NCHAR'] = int(self.objects['DIMENSIONS']['NCHAR'])
        self.__checkCharlabels(nc)

        if self.objects['DIMENSIONS'].has_key('NTAX') :
            self.objects['DIMENSIONS']['NTAX'] = int(self.objects['DIMENSIONS']['NTAX'])
        nitems = len(self.objects['FORMAT']['ITEMS'])
        self.objects['TAXLABELS'] = []
        missing = self.objects['FORMAT']['MISSING']
            
        i = 0
        while i < len(token_list) :
            item = []
            name = token_list[i]
            self.objects['TAXLABELS'].append(name)
            i+=1
            for char in range(nc) :
                if token_list[i] == missing :
                    i+=1 # just skip
                elif token_list[i] == '(' :
                    # TODO: if there are multiple items, they must not have missing values!
                    self.objects['MATRIX'][(name,self.objects["CHARLABELS"][char])]  = (map(float,(token_list[i+1:i+1+nitems])))
                    i = i+2+nitems # skip past ')'
                else :
                    self.objects['MATRIX'][(name,self.objects["CHARLABELS"][char])] = (float(token_list[i]))
                    i+=1
            #self.objects['MATRIX'][name]= item



##############################################################################
## CactusBlock
##############################################################################

class CactusBlock(ContinuousBlock):
    '''CACTUS block.

       Inherits from ContinuousBlock and stores
        program preferences in addition to a simple matrix.
    '''
    def __init__(self, n = 'CACTUS', log = sys.stdout):
        ContinuousBlock.__init__(self, n, log)

        self.objects['PROPERTIES'] = {
            'OVERWRITE_CONTINUOUS' : False,
            'WRITE_BRANCH_LENGTHS' : False,
            'WRITE_SETS'           : True,
            'CHARSET'              : 'CACTUS_CHAR_SET',
            'TREESET'              : 'CACTUS_TREE_SET',
            'PRUNESET'             : 'CACTUS_PRUNE_SET',
            'PRECISION'            : 9}


    def __repr__(self):
        return self.asString()
  
    def asString(self):        
        result = ["BEGIN %s;\n" % self.blockname,]

        result.append(self.writeDimensions())
        result.append(self.writeFormat())
        result.append(self.writeCharlabels())
        result.append(self.writeMatrix())
        
        result.append('PROPERTIES ')
        for p, val in  self.objects['PROPERTIES'].items():
            result.append("%s = %s " % (p,nx_string(val)))
        result.append(';\n')
  
        result.append('END;\n')
        return ''.join(result)

##############################################################################
## TreesBlock
##############################################################################

class TreesBlock(BlockProcessor):
    '''TREES block.

       stores the following data:
        object["TRANSLATE"]: a TRANSLATE table
        object["TREES"]     a dictionary containing the PhyloTree objects
    '''
    def __init__(self, n = 'TREES', log = sys.stdout):
        BlockProcessor.__init__(self, n, log)
        self.objects["TRANSLATE"] = NexusDict()
        self.objects["TREES"] = NexusDict()

    def __repr__(self):
        return self.asString()

    def asString(self, with_translate=False):   
        "String representation of TREES block"
        result = ["BEGIN %s;\n" % self.blockname,]
       # print len(self.objects["TREES"].items())
        
        if with_translate:
            self.make_translate()
            for name, tree in self.objects["TREES"].items():
                tree.relabel_taxa(self.objects["TRANSLATE"])
            result.append("TRANSLATE")
            entries =  self.objects["TRANSLATE"].items()
            for tup in entries[:-1]:
                result.append("%s\t%s," % tup)
            result.append("%s\t%s\n;\n") % entries[-1]
            #result.append(self.objects["TRANSLATE"].__repr__())
            
        # write the trees
        for name, tree in self.objects["TREES"].items():
            result.append("\tTREE %s = %s;" %(nx_string(name), tree.write(True)))
            
        result.append("END;")
        return '\n'.join(result)

    def make_translate(self):
            """Create translate table"""
            tips = sum(map(lambda t:t[1].leaves(),self.objects["TREES"].items()),[])
            taxa = map(lambda n:n.label, tips)
            trans = self.objects["TRANSLATE"]
            for t in taxa : trans[t]=""
            count =1
            for name in trans.keys():
                trans[name]= "%d" % count
                count += 1
                

    def assignObject(self, category, object_name, object_description, format = None, has_star = 0) :
        '''Handle assignments of type: CATEGORY * NAME = (FORMAT) = description.'''
    
        if category == 'TREE' :
            t = newick.create_tree(object_description)
            t.relabel_taxa(self.objects["TRANSLATE"])
            self.objects["TREES"][object_name] = t
            # TODO: record rootedness
        else :
            raise InputError(category, 'Unrecognized object category')
            #print >> self.log, 'Unrecognized object category %s' % category
            
    def doTRANSLATE(self, token_list):
        '''Reads TRANSLATE table.'''
        for i in range(0,len(token_list),3):
            self.objects["TRANSLATE"][ token_list[i] ] =  token_list[i+1]
            #print self.objects["TRANSLATE"]

# end: class TreesBlock

##############################################################################
## ContinuousBlock
##############################################################################
class SetsBlock(BlockProcessor):
    '''SETS block.
       stores sets as list objects.
       to access: self.objects["CHARSET"]["MYSET"]

    '''
    def __init__(self, n = 'SETS', log=sys.stdout):
        BlockProcessor.__init__(self, n, log=log)

    def __repr__(self):
        return self.asString()

    def asString(self):   
        result = ["BEGIN %s;\n" % self.blockname,]
        for cat, obj in self.objects.items():
            for name, s in obj.items():
                if len(s) > 0 :  # don't print empty sets
                    result.append("%s %s (VECTOR) = " % (cat, nx_string(name)))
                    for i in s :
                        result.append("%d " % (i+1))
                    result.append(';\n')
        result.append('END;\n')
        return ''.join(result)
            
 
    def readSet(self, description, format = "VECTOR") :
        "Read set into list.  Maps integer indexes from 1-based to 0-based"
        try :
            if format.upper() == "VECTOR" :
                return map(lambda x : int(x)-1,description)
            else :
                raise InputError(' '.join(description), 'Unrecognized set format: %s' % format)
        except :
            raise InputError(' '.join(description), 'Unrecognized set')

    def assignObject(self, category, object_name, object_description, format = None, has_star = 0) :
        '''Handle assignments of type: CATEGORY * NAME = (FORMAT) = description.'''
        self.addAttribute(category, object_name, self.readSet(object_description, format))



# --------------- Test --------------- #
if __name__ == '__main__' :
    import nexus_grammar, nexus_parser
    reload(nexus_grammar)
    src = open(sys.argv[1]).read()
    
    nx = nexus_parser.NexusParser()
    
    nx.addRecognize('TREES', TreesBlock)
    nx.addRecognize('CONTINUOUS', ContinuousBlock)
    nx.addRecognize('SETS',SetsBlock)

    nexus_grammar.parser.parse(src, processor = nx)
    for k, b in nx.blocks.items():
       print b
      # pass

    #print nx.blocks["TREES"]

    

    
