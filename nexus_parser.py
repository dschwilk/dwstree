#! /usr/bin/env python

# File: nexus_parser.py
# Author: Dylan Schwilk (www.pricklysoft.org)
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

'''
Provides processors for NEXUS files and NEXUS blocks.  The parsers
inherit from SimpleParse DispatchProcessor.
'''

__version__ = "1.3"
__author__= "Dylan Schwilk"

import nexus_grammar, mx.TextTools, sys
from nexus_dict import NexusDict   # these keeps sequence
from simpleparse.dispatchprocessor import *

# constants
ws = mx.TextTools.set(" \t\n\r\'\"")

#Message types
ERROR = "ERROR"
WARNING = "WARNING"
MESSAGE = "MESSAGE"
COMMENT = "COMMENT"

## Exceptions
class NexusError(StandardError):
    """Base class for exceptions in this module."""
    pass

class InputError(NexusError):
    """Exception raised for errors in the input.

    Attributes:
        expression -- input expression in which the error occurred
        message -- explanation of the error
    """

    def __init__(self, expression, message):
        self.expression = expression
        self.message = message


# ------------- NexusFile ----------------#
class NexusParser( DispatchProcessor ):
    '''Dispatch processor for NEXUS grammar to create a NEXUS file.
       Provides file level dispatcher for taglist produced by
       nexus_gramar.py.

       By default, the processor skips all blocks.  To process blocks,
       add block names and classes to the the addRecognize() function
       before parsing.'''
    
    def __init__(self, log = sys.stdout, store_unrecognized = True) :
       
        self.store_unrecognized = store_unrecognized
        self.blocks = NexusDict() # add block objects to this list
        self.comments = []  # list of comments found between blocks
        self.log = log      # where log output should go
        self.__recognize = {} #dict blocks to recognize        

    def __repr__(self) :
        '''Represent NEXUS file as a string for file storage'''
        result = ['#NEXUS\n']
        for key, block in self.blocks.items() :
            if type(block) == type(""):
                result.append(block)  # an unrecognized block is stored as a string
            else :
                result.append(self.blocks[key].asString())
        return '\n'.join(result)    

    def addRecognize(self, name, block_class):
        self.__recognize[name] = block_class

    def logMessage(self, type, message, (start, buffer)):
        if self.log :
            print >> self.log, "Line %d - %s: %s" % (lines(0,  start,  buffer ), type, message)
        

    # ------- simpleparse.dispatchprocessor taglist callbacks -----------#
    def block(self, (tag,start,stop,subtags), buffer ):
        '''Check if block name is in dictionary of blocks to process,
           if so dispatch on all commands in block.'''

        # first skip past leading comments
        map = multiMap(subtags)
        name = (dispatch(self,  map['name'][0] , buffer)).upper()
        if self.__recognize.has_key(name) :
            self.logMessage(MESSAGE, "Processing block '%s'" % name, (start,  buffer))
            self.blocks[name] = self.__recognize[name](log = self.log)  # create new instance of block class
            try :
                dispatchList(self.blocks[name], map['command'], buffer)
            except :
                self.logMessage(ERROR, "Empty or incorrect block '%s'" % name, (start,  buffer))   
        else :
            self.logMessage(MESSAGE, "Skipping block '%s'" % name, (start,  buffer))
            if self.store_unrecognized :
                self.blocks[name] = getString((tag,start,stop,subtags), buffer)



    # ------------- Need these calls to process block names ----------------#
    #               This does duplicate info in NexusBlock class, though  -- cleaner way?
    def word(self, (tag,start,stop,subtags), buffer ):
        '''We need to get the word parts separately because a word can contain comments and strings'''
        return ''.join(dispatchList(self, subtags, buffer))

    def word_part(self, (tag,start,stop,subtags), buffer ):
        '''We need to get the word parts separately'''
        return getString((tag,start,stop,subtags), buffer)

    def string(self, (tag,start,stop,subtags), buffer ):
        '''Need to fix so that it deals with two quotes'''
        return getString((tag,start+1,stop-1,subtags), buffer)

    def name(self, (tag,start,stop,subtags), buffer ):
        return dispatch(self, subtags[0], buffer)

    def comment(self, (tag,start,stop,subtags), buffer ):
        '''Comments outside of a block'''
        c = getString((tag,start,stop,subtags), buffer).strip()
        self.comments.append(c)
        # print to log if output comment
        if c[1] in "!&%/\@" :
            self.logMessage(COMMENT, c, (start,  buffer))
        return ''  # comments stripped from tokens at this level

 
#
# ------------- BlockProcessor ----------------
# Base class for nexus blocks

class BlockProcessor( DispatchProcessor ):
    '''
       Base class dispatch processor for NEXUS block.  Provides basic
       facilities. This class is able to recognize simple assignment
       (object definition) commands where assignment is a word or
       number or.  All other commands must be defined in the derived
       class.  Commands are called according to the following syntax.
       For commands recognized as object assignments::

           addAttribute(category, object_name, object_description, has_star = False)

       is called.  This function should not have to be overridden and simply assingns
       values to the self.objects dictionary.  For other assignments, block classes
       should override::

           assignObject(self, category, object_name, object_description, format = None, has_star = 0)

       For all other commands. doCommand is called which in turn calls the function
       doCOMMANDNAME(list) if such a command exists.


        At the bare minumum, classes derived from NexusBlock must
        implement any commands that the block is expected to contain.  Note that some object
        assignments are not recognized by nexus_parser.py and will end
        up being called as regular commands.'''
        

    def __init__(self, blockname, log = sys.stdout) :
        self.log = log                      # Error output
        self.blockname = blockname
        self.objects = NexusDict() # dict of dicts self.objects[category][object_name] =  (object , has_star)

    def __repr__(self) :
        '''default representation'''
        r = ['BEGIN ;\n' % self.blockname]
        for k in self.objects.keys() :
            r.append( '\n%s %s;\n' % (k, self.objects[k]))
        r.append('END;\n')
        return ''.join(r)

    def asString  (self) :
        return self.__repr__()

    # These are the main public functions used by clients:
    def addAttribute(self, category, object_name, object_description, has_star = 0) :
        '''Default implementation just adds object to
        objects[category].  Should work fine for most blocks.'''     
        if not self.objects.has_key(category) :
            self.objects[category] = NexusDict()
        self.objects[category][object_name] = object_description
        if has_star : self.objects[category]['DEFAULT'] = object_description

    def assignObject(self, category, object_name, object_description, format = None, has_star = 0) :
        '''Handle assignments of type: CATEGORY * NAME = (FORMAT) = description.'''
        self.addAttribute(category, object_name, object_description)
             
    def doCommand(self, command_name, token_list) :
        '''Default implementation calls command do[command_name](token_list).'''
        exec('self.do%s(token_list)' % command_name.upper() )

    def getObject(self, cat, name):
        """Get object, return None if objects does not exist"""
        try :
            return self.objects[cat][name]
        except KeyError :
            return None

   
    # ------- simpleparse.dispatchprocessor taglist callbacks -----------#
    # These methods will not generally be called by clients, only during
    # dispatch processing
    def command( self, (tag,start,stop,subtags), buffer ):
        '''Process command'''
        dispatchList(self, subtags, buffer)

    def simple_assignment(self, (tag,start,stop,subtags), buffer ):
        '''Add object.'''

        map = multiMap(subtags)
        category = dispatch(self, map['category'][0], buffer)
        has_star = 0
        if 'star' in map.keys() : has_star = 1
        for t in map['assignment'] :
            n, d = dispatch(self, t, buffer)
            self.addAttribute(category, n, d, has_star = has_star)

    def complex_assignment(self, (tag,start,stop,subtags), buffer ):
         has_star = 0
         the_map = singleMap(subtags)
         if 'star' in the_map.keys() : has_star = 1
         format = the_map.get('object_format', None)
         if format : format = dispatch(self, format, buffer)
         cat = dispatch(self, the_map['category'], buffer).upper()

         try :
             self.assignObject(cat, dispatch(self, the_map['name'], buffer),
                               dispatch(self, the_map['word_token_list'], buffer),
                               format = format, has_star = has_star)
         except InputError, e :
             self.logMessage(ERROR, "Unrecognized input '%s' - %s" % (e.expression, e.message), (start,  buffer))

         
    def assignment(self, (tag,start,stop,subtags), buffer ):
        '''Return tuple: name, format, description.'''
        the_map = singleMap(subtags)
        return (dispatch(self, the_map['name'], buffer),
                dispatch(self, the_map['object_description'], buffer))

    def object_description(self, (tag,start,stop,subtags), buffer ):
        return dispatch(self, subtags[0], buffer)

    def category(self, (tag,start,stop,subtags), buffer ):
        return dispatch(self, subtags[0], buffer)

    def name(self, (tag,start,stop,subtags), buffer ):
        '''name is just a word so call on word'''
        return dispatch(self, subtags[0], buffer)

    def object_format(self, (tag,start,stop,subtags), buffer ):
        '''format is just a word so call on word'''
        return dispatch(self, subtags[0], buffer)
         
    def other_command( self, (tag,start,stop,subtags), buffer ):
        '''Attempt to call function __do_command with list of results'''
        #ls = []
        cmd = (dispatch(self, subtags[0], buffer)).upper()
        ls = filter(None, dispatchList(self, subtags[1:], buffer))      
        try :
             self.doCommand(cmd, ls)
        except (AttributeError) :
            self.logMessage(ERROR, "Unrecognized command '%s'" % cmd, (start,  buffer))
        except InputError, e  :
            self.logMessage(ERROR,
                            "Command '%s' syntax error at '%s' - %s" % \
                            (cmd, e.expression, e.message), (start,  buffer))

    def token_list(self, (tag,star,stop,subtags), buffer ):
        '''return list of tokens.'''
        return filter(None, dispatchList(self, subtags, buffer))

    def word_token_list(self, (tag,star,stop,subtags), buffer ):
        '''return list of tokens.'''
        return filter(None, dispatchList(self, subtags, buffer))

    def comment(self, (tag,start,stop,subtags), buffer ):
        return ''.join(dispatchList(self, subtags, buffer))
  
    def ignore_comment(self, (tag,start,stop,subtags), buffer ):
        "Non command comments are returned as emtpy strings"
        return ''

    def command_comment(self, (tag,start,stop,subtags), buffer ):
        '''Default implementation returns empty string.
            Output comments are printed to self.log
            Derived classes should override the comment implementation
            when they need to preserve comments'''
        c = getString((tag,start,stop,subtags), buffer).strip()
        #print >> self.log, 'Line %d - %s' % (lines(0,  start,  buffer ),c)
        self.logMessage(COMMENT, c, (start,  buffer))
        return ''  # comments stripped from tokens at this level 


    def word(self, (tag,start,stop,subtags), buffer ):
        '''We need to get the word parts separately because a word can
        contain comments'''
        return ''.join(dispatchList(self, subtags, buffer))

    def word_part(self, (tag,start,stop,subtags), buffer ):
        '''We need to get the word parts separately'''
        return getString((tag,start,stop,subtags), buffer)

    def safepunct(self, (tag,start,stop,subtags), buffer ):
        '''Return punctuation character.'''
        return getString((tag,start,stop,subtags), buffer)            

    def punct_no_equal(self, (tag,start,stop,subtags), buffer ):
        '''Return punctuation character.'''
        return getString((tag,start,stop,subtags), buffer)
    
    def tuple(self, (tag,start,stop,subtags), buffer ):
        '''( token ...)'''
        return filter(None, dispatchList(self, subtags, buffer))  # filter eliminates ignore comments

    # numbers returned as strings.  Client is responsible for converting
    def number(self, (tag,start,stop,subtags), buffer ):
        #return dispatch(self, subtags[0], buffer)
        return (getString((tag,start,stop,subtags), buffer))
    
##    def int(self, (tag,start,stop,subtags), buffer ):
##        return (getString((tag,start,stop,subtags), buffer))
##        #return int(getString((tag,start,stop,subtags), buffer))
##
##    def float(self, (tag,start,stop,subtags), buffer ):
##        return float(getString((tag,start,stop,subtags), buffer))


    def star(self, (tag,start,stop,subtags), buffer ):
        return 1

    def string(self, (tag,start,stop,subtags), buffer ):
        '''Need to fix so that it deals with two quotes'''
        return getString((tag,start+1,stop-1,subtags), buffer)


    def logMessage(self, type, message, (start, buffer)):
        if self.log :
            print >> self.log, "Line %d - %s: %s" % (lines(0,  start,  buffer ), type, message)


# utility functions for writing nexus files

def nx_string(s):
    '''Write NEXUS string from number, string or list'''
    if type(s) == type([]) :
        result = ['(']
        for i in s :
            result.append(' %s' % str(i))
        result.append(')')
        return ''.join(result)
    elif type(s) <> type(""):
        s = str(s)
    
    if mx.TextTools.setfind(s, ws) > -1 :
        # TODO: replace ' and " with '' and ""
        return '"%s"' % s
    return s

