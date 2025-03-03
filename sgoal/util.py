# Useful methods.
# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co  and eleonguz@unal.edu.co
# All rights reserved.
# Licence
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
# in the documentation and/or other materials provided with the distribution.
# Neither the name of the copyright owners, their employers, nor the names of its contributors may be used to endorse or 
# promote products derived from sgoal.this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

import random as rand
from inspect import signature

############### UTILITY FUNCTIONS ################
# Arity of a variation operator f. Used for computing number of required parents
def arity(f): return len(signature(f).parameters)

# Generates a boolean value according to probability p ( True with probability p, False otherwise )
def randbool(p=0.5):
  return (rand.random() < p)

# A permutation of n elements
def permutation(n):
  x = [i for i in range(0,n)]
  rand.shuffle(x)
  return x

# Normalizes a vector of weights
def normalize(weight):
  c = 0
  for w in weight: 
    c += w
  return [ w/c for w in weight ]

# Binary search of value in sorted array
def index(x, A, ini=0, fin=-1):
  if(fin==-1): fin=len(A)
  while(ini+1<fin):
    m = (fin+ini)//2
    if(A[m]==x): return m
    if(A[m]<x): ini = m+1
    else: fin = m
  if(ini==fin or A[ini]!=x): return -1
  return ini

# Power law number generation
def powerlaw(a=-2.0):
  if a==-2.0: return 1/(1-rand.random())
  else: return (1-rand.random())**(1/(a+1))