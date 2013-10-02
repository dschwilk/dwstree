#! /usr/bin/env python

# File: nexus_grammar.py
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
:Author: Dylan W. Schwilk
:Date:   2005/03/09

EBNF grammar and parser for NEXUS file language.  For use with
SimpleParse. This grammar does not try to recognize all constructs but
does more parsing than a simple tokenizer.  It will recognize all
constructs at least to the 'command' level.  But other than simple
object assignments, such commands will simple be parsed as a command
name and a list of tokens.  The NexusBlock object is responsible for
dealing with such commands.

Features:

1. Recognizes nested comments
2. returns adjacent tokens, strings and comments as single word for processor to concatinate, clean
3. Processor class can decide whether to recognize and save comments, whether to strip
   command comments from tokens, etc


Limitations:

1. Does fairly minimal structural parsing/recognition.  Simple objects
   are recognized, other  commands a parsed as token lists and left to
   the processor.  
   
2. Will not recognize mixed-case 'BEGIN' or END' commands; only easy
   way in simple parse would be to hard-code alternatives
   
3. There are certainly BUGS!  Need to test.
"""

from simpleparse.parser import Parser
import pprint
from simpleparse.common import numbers


dec = r'''
nexus_file           := whitespacechar*, NEXUS ,  whitespacechar+, (comment / block)+ 
<NEXUS>              := ! , '#NEXUS' 

# --------------- Blocks --------------#
block           := sep, BEGIN,  sep_or_comment, name,  sep_or_comment, ';' ,
                   sep_or_comment, block_body,  sep_or_comment, END,  sep_or_comment, ';'
command         := simple_assignment / complex_assignment / other_command
other_command   := command_name, (sep_or_comment, token)+ , sep_or_comment , ';'
>command_name<  := word
>BEGIN<         := 'BEGIN' / 'begin'
>END<           := 'END' / 'end' / 'ENDBLOCK' / 'endblock'
>block_body<    := (sep_or_comment, (comment  / command))*

# --------------- Tokens --------------#
name             :=  word
>token<          :=  nx_number / word / safepunct / comment  # This order works!
>word_token<     :=  word / safepunct / comment   # this is to avoid numbers being recognized when not wanted
token_list       :=  (sep_or_comment, token)+
word_token_list  :=  (sep_or_comment, word_token)+
>sep<            := whitespacechar*             # should I allow comments to act as token separators?
>sep_or_comment< := ( whitespacechar / standalone_comment)*
<whitespacechar> := [ \t\n\r]
<ws>             := whitespacechar+
safepunct          := punct_no_equal / [=]
punct_no_equal     := [-()\\{}/,:*+<>]              # so we can recognize assignments
punct              := safepunct / [][;'"]
>nx_number<        := number #, ?[wordchars]           # nexus words can begin with digits, so number must not end with wordchars
<wordchars>      := -(punct/whitespacechar)
word_part        := wordchars+
word             := string / (word_part / comment / string)+   # comments and strings do not break words!
string           := dq_string / sq_string
<dq_string>      := '"', -["]* , '"'
<sq_string>      :=  ['], ?-['] , (-['] / '\'\'' )* , ['], ?-[']

# ---------------  Objects ----------------
# a few object types are explicitly recognized
# all others can be handled by dispatcher
# -----------------------------------------
simple_assignment  := category, sep_or_comment, star?, (sep_or_comment, assignment)+ , sep_or_comment , ';'
complex_assignment := category, (sep_or_comment, star)?, sep_or_comment, name, sep_or_comment,
                       ( [(] , object_format,  ')' , sep_or_comment)?,  '=' , sep_or_comment, word_token_list, sep_or_comment , ';'
assignment         := name, sep_or_comment,  '=' , sep_or_comment, object_description+
object_description :=  nx_number / word / tuple
star               := '*'
category           := word
object_format      := word
tuple              := '(' , sep_or_comment, ( (nx_number / word), sep_or_comment)+ , ')'

# --------------- Comments   ----------------
comment               := sep, (command_comment / ignore_comment)
>standalone_comment< := comment, whitespacechar            # Needed so comments adjacent to tokens are returns as token
command_comment       := '[' , [!&%/\@] , comment_string , ']'
ignore_comment        := '[', comment_string , ']'
<comment_string>      := (-[][]+ / nested_comment)*
<nested_comment>      := ('[' , comment_string, ']')
'''


# --------------- The parser ----------------
parser = Parser(dec,'nexus_file')

# --------------- Test function ----------------
if __name__ == "__main__":
    import sys
    #src = sys.stdin.read()
    src = sys.stdin.read()
    taglist = ( parser.parse( src))
    pprint.pprint(taglist)

    print 'Nexus file has %d blocks or out-of-block comments.' % len(taglist[1])

