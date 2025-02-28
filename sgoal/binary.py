# Binary (BitArray) search space definitions
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
from sgoal.core import SPACE
from sgoal.core import PROBLEM
from sgoal.core import randbool
from sgoal.core import VRSGoal
from sgoal.core import variation
from sgoal.core import simplexover
from sgoal.core import transposition
from sgoal.es import Rule1_5_T
from sgoal.ga import SSGA_T
from sgoal.ga import GGA_T
from sgoal.chavela import CHAVELA_T

# Fixed Length Binary Space
# D : Length of the BitArray (Binary string)
def Binary(D): 
  space = SPACE( lambda: [1 if randbool() else 0 for i in range(D)] )
  space['D'] = D
  return space

############### VARIATION OPERATIONS ################
# Flips the kth bit of BitArray x (creates a new one - hard copy)
def flip(x, k):
  y = x.copy()
  y[k] = 1 if y[k]==0 else 0
  return y

# Generates the complement BitArray by flipping each bit (creates a new one - hard copy)
def complement(x): return [1 if v==0 else 0 for v in x]

# Single bit mutation: Flips a single bit chosen in a random fashion (creates a new one - hard copy)
def singlebitmutation(x):
  return flip(x, rand.randint(0,len(x)-1))

# Multiple bits mutation: Flips the set of bits in the indices array k (creates a new one - hard copy)
def multiflip(x, k):
  if(len(k)==0): k=[i for i in range(len(x))]
  y = x.copy()
  for i in k:
    y[i] = 1 if y[i]==0 else 0
  return y

# Bit mutation. Flips each bit with probability p (creates a new one - hard copy)
def bitmutationprob(x, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): 
      y[i] = 1 - y[i]
  return y

# Bit mutation. Flips a bit with probability 1/|x| (creates a new one - hard copy)
def bitmutation(x): return bitmutationprob(x, 1.0/len(x))

##################### SGOALs ###########################
# Classical Hill Climbing Algorithm with neutral mutation for BitArray problems. Uses bitmutation as variation operator
# problem: Problem to solve
def HC(problem): 
  if( 'variation' not in problem ): problem['variation'] = bitmutation 
  return VRSGoal(problem)

# The HC algorithm suggested by Richard Palmer, that Forrest and
# Mitchell named as "random mutation hill-climbing" (RMHC), see
# M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?”
# Santa Fe Institute, Working Papers, 01 1993.
# problem: Problem to solve
def RMHC(problem): 
  if('variation' not in problem): problem['variation'] = singlebitmutation 
  return VRSGoal(problem)

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, for BitArray
# problem: Problem to solve
def setprob(problem):
  problem['variation'] = lambda x, fx: variation(x, fx, lambda y: bitmutationprob(y, problem['parameter']), problem)
  
# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
def Rule1_5(problem):
  D = problem['D']
  if( 'parameter' not in problem ): problem['parameter'] = 1/D
  if( 'variation' not in problem ): 
    problem['setparameter'] = lambda : setprob(problem)
    setprob(problem)
  if( 'G' not in problem ): problem['G'] = D
  return Rule1_5_T(problem)

############### Generational Genetic Algorithm - SSGA ################
def bmutation(sgoal):
  if('mutation' not in sgoal): sgoal['mutation'] = bitmutation
  return sgoal

def GGA(problem):
  return GGA_T(bmutation(problem))

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
def SSGA(problem):
  return SSGA_T(bmutation(problem))

# Standard CHAVELA for Binary problems. Uses bitmutation, simplexover, and transposition as operators
def CHAVELA(problem):
  if( 'operators' not in problem ): problem['operators'] = [bitmutation, simplexover, transposition]
  return CHAVELA_T(problem) 

# Standard CHAVELA1 for Binary problems
def CHAVELA1(problem):
  if( 'operators' not in problem ): problem['operators'] = [bitmutation, transposition]
  return CHAVELA_T(problem) 

##################### TEST FUNCTIONS #####################
# Computing the MaxOnes function (sum of bits) from the start bit upto end-1 bit. 
# Sums bits up to the end of the BitArray if end is set to -1
def maxones(x, start=0, end=-1):
  if(end<0): 
    end = len(x)
  s = 0
  for i in range(start,end): 
    s += x[i]
  return s

# Goldberg's 3-Deceptive function, as defined in 
# D. Goldberg, B. Korb, and K. Deb, “Messy genetic algorithms: 
# motivation, analysis, and first results,” 
# Complex Systems, vol. 3, pp. 493–530, 1989.
def deceptive(x, start=0, end=-1):
  if(end<0): 
    end = len(x)
  switcher = {
    0: 28,
    1: 26,
    2: 22,
    3: 0,
    4: 14,
    5: 0,
    6: 0,
    7: 30
  }
  size=3
  c = 0
  i = start
  while i<end:
    d=0
    b=1
    for k in range(i,i+size):
      d += b if x[k] else 0
      b *= 2
    c += switcher.get(d, 0)
    i += size  
  return c

# Goldberg's Boundedly-Deceptive function, as defined in 
# D. Goldberg, B. Korb, and K. Deb, “Messy genetic algorithms: 
# motivation, analysis, and first results,” 
# Complex Systems, vol. 3, pp. 493–530, 1989.
def generic_boundedly(x, size, start=0, end=-1):
  if(end<0): 
    end = len(x)
  c = 0
  i = start
  while i<end:
    u = 0
    for k in range(i,i+size):
      u += 1 if x[k] else 0
    c += u if u==size else size-1-u
    i += size
  return c

def boundedly(x, start=0, end=-1): return generic_boundedly(x,4,start,end)

# Forrest's Royal Road function as defined in
# M. Mitchell, S. Forrest, and J. Holland, “The royalroad 
# for genetic algorithms: Fitness landscapes and ga performance,” 
# in Toward a Practice of Autonomous Systems: 
# Proceedings of the First European Conference on Artificial Life, 11 1992.
def generic_royalroad(x, size, start=0, end=-1):
  if(end<0): 
    end = len(x)
  c = 0
  i = start
  while i<end:
    start += size
    while i<start and x[i]:
      i+=1
    c += size if i==start else 0
    i = start
  return c

def royalroad8(x,start=0,end=-1): return generic_royalroad(x,8,start,end)

def royalroad16(x,start=0,end=-1): return generic_royalroad(x,16,start,end)

# The mixed function as defined in 
# J. Gomez and E. Leon, "Gabo: Gene Analysis Bitstring Optimization," 
# 2022 IEEE Congress on Evolutionary Computation (CEC), Padua, Italy, 2022, pp. 1-8, doi: 10.1109/CEC55065.2022.9870237.
# A combination of all the previous bit functions, each block of 20 bits is defined as follow
# 0..4 : maxones
# 5..7 : deceptive3
# 8..11 : boundedly
# 12..19 : royalroad8
def mixed(x, start=0, end=-1):
  if(end==-1): 
    end = len(x)
  f = 0
  while(start<end):
    f += maxones(x,start,start+5) + deceptive(x,start+5,start+8) + boundedly(x,start+8,start+12) + royalroad8(x,start+12,start+20) 
    start += 20
  return f 

# The mixed function 2: Defined in a similar fashion to mixed but using RR2 (royalroad 16  bits)
# See, J. Gomez and E. Leon, "Gabo: Gene Analysis Bitstring Optimization," 
# 2022 IEEE Congress on Evolutionary Computation (CEC), Padua, Italy, 2022, pp. 1-8, doi: 10.1109/CEC55065.2022.9870237.
# A combination of all the previous bit functions, each block of 20 bits is defined as follow
# 0..9 : maxones
# 10..15 : deceptive3
# 16..23 : boundedly
# 23..39 : royalroad16
def mixed2(x, start=0, end=-1):
  if(end==-1): 
    end = len(x)
  f = 0
  while(start<end):
    f += maxones(x,start,start+10) + deceptive(x,start+10,start+16) + boundedly(x,start+16,start+24) + royalroad16(x,start+24,start+40) 
    start += 40
  return f 

##################### TEST PROBLEMS ####################
def TestProblem(f, D, EVALS, TRACE=False):
  space = Binary(D)
  space['optimum'] = D
  if(f=='MaxOnes'): f = maxones
  elif(f=='GD3'): 
    f = deceptive
    space['optimum'] = 10*D
  elif(f=='GD4'): f = boundedly
  elif(f=='RR1'): f = royalroad8
  elif(f=='RR2'): f = royalroad16
  elif(f=='Mixed'):
    f = mixed
    space['optimum'] = 47*D/20
  elif(f=='Mixed2'):
    f = mixed2
    space['optimum'] = 47*D/20
  else: f = maxones
  return PROBLEM('max', f, space, EVALS, TRACE)
