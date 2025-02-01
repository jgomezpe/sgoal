# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Binary search space definitions

import random as rand
from sgoal import Space, transposition
from sgoal import SGoal
from sgoal import tournament
from sgoal import simplexover
from sgoal import randbool
from hc import HC
from ga import GGA
from ga import SSGA
from chavela import CHAVELA
from rule_1_5th import Rule_1_5th

# BitArray Space
class BitSpace(Space):
  def __init__(self, D):
    self.D = D

  def getone(self):
    return [1 if randbool() else 0 for i in range(self.D)]
  
############### USEFUL FUNCTIONS ################

# For all logic operator tested on an array of boolean values
def for_all(a):
  for v in a:
    if(not v):
      return False
  return True

############### VARIATION OPERATIONS ################
# Flips the kth bit of a genome (creates a new genome)
def flip(x, k):
  y = x.copy()
  y[k] = 1 if y[k]==0 else 0
  return y

# Generates the complement bitstring (creates a new genome)
def complement(x): return [1 if v==0 else 0 for v in x]

# Single bit mutation (creates a copy with a bit randomly flipped)
def singlebitmutation(x):
  return flip(x, rand.randint(0,len(x)-1))

# Multiple bit mutation
def multiflip(x, k):
  y = x.copy()
  for i in k:
    y[i] = 1 if y[i]==0 else 0
  return y

# Bit mutation. Flips a bit with probability p
def bitmutationprob(x, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): 
      y[i] = 1 - y[i]
  return y

# Bit mutation. Flips a bit with probability 1/|x|
def bitmutation(x): return bitmutationprob(x, 1.0/len(x))


#################### CLASSIC LITERATURE ALGORITHMS #####################
# The HC algorithm suggested by Richard Palmer, that Forrest and
# Mitchell named as "random mutation hill-climbing" (RMHC), see
# M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?”
# Santa Fe Institute, Working Papers, 01 1993.
# problem: Problem to be solved
def RMHC(problem): return HC(problem, singlebitmutation)


# A Generalization of The Global Search Algorithm for GA-easy functions 
# propossed by Das and Whitley 1991, which tries only order 1-schemas, see
# R. Das and L. D. Whitley, “The only challenging problems are deceptive: 
# Global search by solving order-1 hyperplanes,
# in Proceedings of the 4th International Conference on Genetic Algorithms, San Diego, CA, USA,
# July 1991 (R. K. Belew and L. B. Booker, eds.), pp. 166– 173, Morgan Kaufmann, 1991.
# Our generalization allows to check a good 1-schema after some CHECK evaluations
# instead of checing it at the end of the allowed number of fitness evaluations
# f: Function to be optimized
class GGS1(SGoal):
  def __init__(self, problem, CHECK):
    SGoal.__init__(self, problem)
    self.S1 = []
    self.fS1 = []
    self.N = 1
    self.check = CHECK
    self.delta = 1

  def evalone(self, x):
    fx = SGoal.evalone(self, x)
    self.S1.append(x)
    self.fS1.append(fx)
    return fx

  def next(self, P, fP):
    if(self.check == -1): 
      self.check = self.evals
    if((self.count+self.delta)%self.check != 0):
      y = self.space.get(1)
      return y, self.evalone(y)

# Computes schemata information
    M = len(self.S1)
    D = self.space.D
    C = [[0 for k in range(D)], [0 for k in range(D)]]
    fH = [[0 for k in range(D)], [0 for k in range(D)]]
    for k in range(D):
      for i in range(M):
        fH[self.S1[i][k]][k] += self.fS1[i]
        C[self.S1[i][k]][k] += 1
    # Generates a candidate solution with the best genes
    y = []
    for k in range(D):
      if( self.minimize ):
        y.append( 1 if(fH[1][k]/C[1][k] < fH[0][k]/C[0][k]) else 0 )
      else:
        y.append( 1 if(fH[1][k]/C[1][k] > fH[0][k]/C[0][k]) else 0 )
    fy = SGoal.evalone(self, y)
    return y, fy

def GS1(problem): return GGS1(problem, -1)

def DGS1(problem): return GGS1(problem, problem['space'].D)


# A Generalization of The Global Search Algorithm with complement propossed by 
# G. Venturini 1995
# Applies the GS1 and compares the obtained solution with its complement and the
# best candidate solution of the S1 set (same as the best solution found), see
# G. Venturini, “Ga consistently deceptive functions are not challenging problems”
# in First International Conference on Genetic Algorithms in Engineering Systems: 
# Innovations and Applications, pp. 357–364, 1995.
# Our generalization checks the complement after some fitness evaluations
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Initial point, a random point. It is required for defining the dimension of the space
# fx: f value at point x (if provided)
class GGSC1(GGS1):
  def __init__(self, problem, CHECK):
    GGS1.__init__(self, problem, CHECK)
    self.delta = 2

  def next(self, P, fP):
    y, fy = GGS1.next(self, P, fP)
    if((self.count+1)%self.check == 0):
      y = complement(self.result['x'])
      fy = SGoal.evalone(self, y)
    return y, fy

def GSC1(problem): return GGSC1(problem, -1)

def DGSC1(problem): return GGSC1(problem, problem['space'].D)

def bitGAconfig(D):
  return {'selection': tournament, 'xover': simplexover, 'xr':0.7, 'mutation':bitmutation, 'N':D//2}
  
def bitGGA(problem):
  return GGA(problem, bitGAconfig(problem['space'].D))

def bitSSGA(problem):
  return SSGA(problem, bitGAconfig(problem['space'].D))

def bitCHAVELA(problem):
  return CHAVELA(problem, {'operators':[bitmutation, simplexover, transposition], 'N':problem['space'].D//2}) 

def bitR1_5(problem):
  D = problem['space'].D
  return Rule_1_5th(problem, {'mr':1/D, 'mutation': bitmutationprob, 'G':D, 'a':0.9}) 

##################### TEST FUNCTIONS #####################
# MaxOnes function
def maxones(x, start=0, end=-1):
  if(end<0): 
    end = len(x)
  s = 0
  for i in range(start,end): 
    s += x[i]
  return s

def twopowermaxones(x, start=0, end=-1):
  return abs(2**maxones(x, start, end)-16)

# Goldberg's 3-Deceptive function
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

# Goldberg's Boundedly-Deceptive function
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

# Forrest's Royal Road function
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

  # A combination of all the previous bit functions, each block of 20 bits is defined as follow
# 0..4 : twopowermaxones
# 5..7 : deceptive3
# 8..11 : boundedly
# 12..19 : royalroad8
# 20..24 : twopowermaxones
def mixed2(x, start=0, end=-1):
  if(end==-1): 
    end = len(x)
  f = 0
  while(start<end):
    f += twopowermaxones(x,start,start+5) + deceptive(x,start+5,start+8) + boundedly(x,start+8,start+12) + royalroad8(x,start+12,start+20) 
    start += 20
  return f 

##################### TEST PROBLEMS ####################
def binaryproblem(f, D):
  if(f=='maxones'):
    return {'f':maxones, 'space': BitSpace(D), 'optimum':D, 'type':'max'}
  if(f=='deceptive'):
    return {'f':deceptive, 'space': BitSpace(D), 'optimum':10*D, 'type':'max'}
  if(f=='boundedly'):
    return {'f':boundedly, 'space': BitSpace(D), 'optimum':D, 'type':'max'}
  if(f=='royalroad8'):
    return {'f':royalroad8, 'space': BitSpace(D), 'optimum':D, 'type':'max'}
  if(f=='mixed'):
    return {'f':mixed, 'space': BitSpace(D), 'optimum':47*D/20, 'type':'max'}
  return {'f':maxones, 'space': BitSpace(D), 'optimum':D, 'type':'max'}