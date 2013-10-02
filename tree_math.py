#! /usr/bin/env python

# File: tree_math.py
# Author: Dylan Schwilk (www.pricklysoft.org)
# $Date: 2005/03/09 00:53:50 $
# $Revision: 1.2 $
# $Source: /home/schwilk/code/python/cactus-pie/nexus/RCS/nexus_math.py $
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


"""Module for math functions

   Functions:
     - sample_rep
     - sample_wr
     - midpoint
     - median
     - average
"""

__version__ = "$Revision: 1.2 $"
__author__  = '''Dylan Schwilk (www.pricklysoft.org)'''

import random
import itertools
from copy import copy

def sample_rep(population, k=None):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    if not k : k = n
    _random, _int = random.random, int  # speed hack 
    result = [None] * k
    for i in xrange(k):
        j = _int(_random() * n)
        result[i] = population[j]
    return result

def sample_wr(population, k):
    "Chooses k random elements (with replacement) from a population"
    n = len(population)
    _random, _int = random.random, int  # speed hack 
    return [_int(_random() * n) for i in itertools.repeat(None, k)]




##def sample_rep(l, N=None):
##    "return list of same lengthof length N sampling with replacement"
##    if not N : N = len(l)
##    result = []
##    for i in range(N): result.append(random.choice(l))
##    return result
    
def midpoint(l):
    return float(min(l)+max(l))/2.0
    
def median(l):
    l2 = copy(l)
    l2.sort()
    N = len(l2)
    if N % 2 == 0 :
        return (l2[N/2] + l2[N/2+1]) / 2.0
    else :
        return l2[N/2]

def average(l):
    return sum(l) / float(len(l))
    
def weighted_average(l, weights):
    " Return weighted avereage of l"
    wa =  [i*j for i, j in zip( l, weights)]  
    return  sum(wa) / sum(weights)


class Quadratic:
    """Holds a, b and c parameters to a quadratic function."""
    def __init__(self, a = 0, b= 0, c=0) :
        self.a = float(a)
        self.b = float(b)
        self.c = float(c)

    def __repr__(self):
        return "%fx^2 + %fx + %f = 0" % (self.a, self.b, self.c)

    def solve(self) :
        """Returns two roots of the quadratic function"""
        from math import sqrt
        root = sqrt(self.b*self.b - 4.0*self.a*self.c)
        return  (-self.b + root) / (2.0*self.a) ,  (-self.b - root) / (2.0*self.a)

    def dsolve(self):
        """returns the solution to the derivative of the quadratic"""
        return self.b / (-2.0*self.a)

class Range:
    def __init__(self, seq=None):
        if seq:
            self.min, self.max = min(seq), max(seq)

    def __repr__(self):
        return "(%s,%s)" % (str(self.min) , str(self.max))

    def add(self, num) :
        if seq:
           self.max =  max(num, self.max) 
           self.min = min(num, self.min) 
        else :
            self.min = num
            self.max = num
           

def intersect(a,b) :
    return Range(max(a.min,b.min), min(a.max,b.max))


def median_of_3(a,b,c) :
    """Return median of three numbers"""
    if max(a,b,c) == a :
        return max(b,c)
    elif min(a,b,c) == a : 
        return min(b,c)
    else :
        return a

def get_range(r1,r2) :
    """ """
    x = median_of_3(r1.min, r1.max, r2.min)
    y = median_of_3(r2.max, r2.min, r1.max)
    return Range( (min(x,y), max(x,y)) )



def main():
    b = Range((52,170.1))
    a = Range((50.1, 151))
    print a, b
    print get_range(a,b)


    q = Quadratic(3,7,1)
    print q
    print q.solve()
    print q.dsolve()

    x = q.solve()[1]
    assert  q.a*(x*x) + q.b*x + q.c == 0.0

    
# Main Test function               
if __name__ == '__main__':

    main()
    
  
