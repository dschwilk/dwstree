dwstree
=======

dwstree provides a python library for reading and writing newick trees and NEXUS files, manipulating phylogenetic trees, and conducting analysis of phylogenies and traits.

Command-line programs available
-------------------------------

All of these modules provide functions/methods and also provide a command-line
interface. All of these programs read data in newick (*.new) and/or NEXUS file
format. For instructions and a list of options, use the '-h' option when
calling the program. For example: type 'python dorder.py -h'.

    - icontrasts.py:     Calculates independent contrasts.
    - dorder.py:         Provides the Divergence Order Test (DOT) and
                         Synchronized Changes Test (SvS).
    - branch_lengths.py: Assign branch lengths to a phylogeny
                         according to several different possible
                         algorithms.
    - cladelabel.py      label clades in a phylogeny as defined by a taxa
                         list in an auxilary file.


Library functions
-----------------

Use the Python help command or browse the code to understand the phylogenetic
tree functions

Installation
------------

> python setup.py install
