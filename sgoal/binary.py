# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Binary search space definitions

import random as rand
from sgoal.util import randbool
from sgoal.core import pick
from sgoal.core import evaluate
from sgoal.core import MAXIMIZE
from sgoal.hc import HC

##################### BIT STRING SOLUTION AND POPULATION #####################
# A bitstring of length L
def bitstring(L):
  return [1 if randbool() else 0 for i in range(L)]

# A population of N bitstrings each one of length L
def bitstring_population(N, L):
  return [bitstring(L) for i in range(N)]
  
############### USEFUL FUNCTIONS ################
# Flips the kth bit of a genome (creates a new genome)
def flip(x, k):
  y = x.copy()
  y[k] = 1 if y[k]==0 else 0
  return y

# For all logic operator tested on an array of boolean values
def for_all(a):
  for v in a:
    if(not v):
      return False
  return True

############### VARIATION OPERATIONS ################
# Generates the complement bitstring (creates a new genome)
def complement(x): return [1-v for v in x]

# Single bit mutation (creates a copy with a bit randomly flipped)
def single_bit_mutation(x):
  return flip(x, rand.randint(0,len(x)-1))

# Bit mutation. Flips a bit with 1/|x| probability
def bit_mutation_probability(x, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): y[i] = 1 - y[i] 
  return y

# Bit mutation with probability. Flips a bit with the given probability
def bit_mutation(x): return bit_mutation_probability(x, 1.0/len(x))

##################### TEST FUNCTIONS #####################
# MaxOnes function
def maxones(x, start=0, end=-1):
  if(end<0): end = len(x)
  s = 0
  for i in range(start,end): s += x[i]
  return s

# Goldberg's 3-Deceptive function
def deceptive(x, start=0, end=-1):
  if(end<0): end = len(x)
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

# Goldberg's Boundedly-Deceptive function
def generic_boundedly(x, size, start=0, end=-1):
  if(end<0): end = len(x)
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

# Forrest's Royal Road function
def generic_royalroad(x, size, start=0, end=-1):
  if(end<0): end = len(x)
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

# A combination of all the previous bit functions, each block of 20 bits is defined as follow
# 0..4 : maxones
# 5..7 : deceptive3
# 8..11 : boundedly
# 12..19 : royalroad8
def mixed(x, start=0, end=-1):
  if(end==-1): end = len(x)
  f = 0
  while(start<end):
    f += maxones(x,start,start+5) + deceptive(x,start+5,start+8) + boundedly(x,start+8,start+12) + royalroad8(x,start+12,start+20)
    start += 20
  return f 
 
#################### CLASSIC LITERATURE ALGORITHMS #####################
# The HC algorithm suggested by Richard Palmer, that Forrest and
# Mitchell named as "random mutation hill-climbing" (RMHC), see 
# M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?” 
# Santa Fe Institute, Working Papers, 01 1993.
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Initial point
# fx: f value at point x
def RMHC( f, evals, x, fx=None): return HC(f, evals, single_bit_mutation, x, fx )

# The Global Search Algorithm for GA-easy functions propossed by Das and Whitley 1991
# Tries only order 1-schemas, see
# R. Das and L. D. Whitley, “The only challenging problems are deceptive: Global search by solving order-1 hyperplanes,
# in Proceedings of the 4th International Conference on Genetic Algorithms, San Diego, CA, USA,
# July 1991 (R. K. Belew and L. B. Booker, eds.), pp. 166– 173, Morgan Kaufmann, 1991.
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Initial point, a random point. It is required for defining the dimension of the space 
# fx: f value at point x (if provided)
def GS1( f, evals, x, fx=None ):
  D = len(x) # Space dimension
  if( evals>0 and not fx):
    fx = evaluate(f, x)
    evals -= 1
  S1 = [x]
  fS1 = [fx]
  for i in range(evals-1):
    x = bitstring(D)
    fx = evaluate(f,x)
    S1.append(x)
    fS1.append(fx)
  # Computes schemata information
  M = len(S1)
  C = [[0 for k in range(D)], [0 for k in range(D)]]
  fH = [[0 for k in range(D)], [0 for k in range(D)]]
  for k in range(D):
    for i in range(M):
      fH[S1[i][k]][k] += fS1[i]
      C[S1[i][k]][k] += 1
  # Generates a candidate solution with the best genes
  y = []
  for k in range(D):
    if( MAXIMIZE ):
      y.append( 1 if(fH[1][k]/C[1][k] > fH[0][k]/C[0][k]) else 0 )
    else:
      y.append( 1 if(fH[1][k]/C[1][k] < fH[0][k]/C[0][k]) else 0 )      
  return y, evaluate(f,y), S1, fS1

# The Global Search Algorithm with complement propossed by G. Venturini 1995
# Applies GS1 and compares the obtained solution with its complement and the 
# best candidate solution of the S1 set (same as the best solution found), see 
# G. Venturini, “Ga consistently deceptive functions are not challenging problems”
# in First International Conference on Genetic Algorithms in Engineering Systems: Innovations and Applications, pp. 357–364, 1995.
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Initial point, a random point. It is required for defining the dimension of the space 
# fx: f value at point x (if provided)
def GSC1( f, evals, x, fx=None ):
  x, fx, S1, fS1 = GS1(f, evals-1, x, fx)
  xc = complement(x)
  fxc = evaluate(f,xc)
  x, xc, fx, fxc = pick(x, xc, fx, fxc)
  for i in range(len(S1)):
    x, y, fx, fy = pick(x, S1[i], fx, fS1[i])
  return x, fx