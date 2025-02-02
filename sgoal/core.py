# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Stochastic Global Optimization Algorthms (SGoals) library core functions.
# A formal description is provided by Jonatan Gomez in
# "Stochastic global optimization algorithms: A systematic formal approach",
# Information Sciences, Volume 472, 2019, Pages 53-76, ISSN 0020-0255,
# https://doi.org/10.1016/j.ins.2018.09.021.
# (https://www.sciencedirect.com/science/article/pii/S0020025517305248)
import random as rand
from inspect import signature

############### UTILITY FUNCTIONS ################
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

def best(quality, MINIMIZE=True):
  k = 0
  if(MINIMIZE):
    for i in range(1,len(quality)):
      if(quality[i] <= quality[k]): 
        k = i
  else:
    for i in range(1,len(quality)):
      if(quality[i] >= quality[k]): 
        k = i
  return k

# Arity of a function f
def arity(f): return len(signature(f).parameters)
  
############### SELECTION MECHANISMS ################
# Uniform selection. Picks N elements (indices) with the same probability
def uniform(quality, N, MINIMIZE=True):
  return [rand.randint(0,len(quality)-1) for i in range(N)]

# Tournament selection. Selects n individuals.
# For each individual to be selected this method picks 4 individuals and returns the best of those
def tournament1(quality, MINIMIZE=True):
  m = 4 #Tournament's size
  candidate = uniform(quality, m, MINIMIZE)
  x = 0
  for k in range(1,m):
    if (MINIMIZE and quality[candidate[x]] >= quality[candidate[k]]) or (not MINIMIZE and quality[candidate[x]] <= quality[candidate[k]]):
      x = k
  return candidate[x]
  
def tournament(quality, N, MINIMIZE=True):
  return [tournament1(quality,MINIMIZE) for i in range(N)]

# Weighted selection: p_i is the probability of selecting element i
def weighted(p):
  y = rand.random()
  k=0
  while k<len(p) and y>=p[k]:
    y -= p[k]
    k+=1
  return k
  
# Roulette wheel selection. Selects n individuals.
def roulette(quality, N, MINIMIZE=True):
  weight = [-q if MINIMIZE else q for q in quality]
  p = normalize(weight)
  return [weighted(p) for i in range(N)]

############### SEARCH SPACE ################
class Space:
  def getone(self):
    return []

  def get(self, N):
    if(N==1):
      return self.getone()
    return [self.getone() for i in range(N)]

############### VARIATION OPERATIONS ################
# Simple crossover.
def simplexover( x1, x2 ):
  p = rand.randint(1,len(x1)-1)
  y1 = x1[0:p] + x2[p:len(x2)]
  y2 = x2[0:p] + x1[p:len(x1)]
  return y1, y2

# Transposition
def transposition(x):
  x = x.copy()
  start = rand.randint(0,len(x)-1)
  end = rand.randint(0,len(x)-1)
  if start>end: start, end = end, start
  while start<end:
    x[start], x[end] = x[end], x[start]
    start += 1
    end -= 1
  return x

############### STOCHASTIC GLOBAL OPTIMIZATION ALGORITHM ################
class SGoal:
  def __init__(self, problem):
    self.f = problem['f']
    self.space = problem['space']
    self.minimize = problem['type'].lower()=='min'
    self.optimum = problem['optimum'] if 'optimum' in problem else None
    self.count = 0
    self.result = {}
    self.TRACE = False
    self.trace = []
    self.poptrace = []
    self.besttrace = []
    self.N = 1

  def pick(self, x, fx, y, fy):
    if((self.minimize and fy <= fx) or (not self.minimize and fy >= fx)):
        return y, fy, x, fx
    return x, fx, y, fy

  def evalone(self, x):
    self.count+=1
    fx = self.f(x)
    if('x' not in self.result):
      self.result['x'] = x
      self.result['f'] = fx
      self.result['evals'] = 1
    else:
      self.result['x'], self.result['f'], b, fb = self.pick(x, fx, self.result['x'], self.result['f'])
      if(x == self.result['x'] and fx != fb):
        self.result['evals'] = self.count
    if(self.TRACE):
      self.trace.append(fx)
      self.besttrace.append(self.result['f'])
    return fx

  def evalpop(self, P):
    i=0
    fP = []
    while(self.caneval() and i<len(P)):
      fP.append(self.evalone(P[i]))
      i+=1
    return fP

  def eval(self, x):
    if(self.N==1):
      return self.evalone(x)
    return self.evalpop(x)

  def caneval(self, n=1):
    return self.count+n <= self.evals
  
  def optimumfound(self):
    return self.result['f']==self.optimum
   
  def stop(self):
    return not self.caneval() or self.optimumfound()

  def next(self, P, fP):
    return P, fP

  def tracing(self, P, fP):
    if(self.TRACE):
      if(self.N > 1):
        if(self.minimize): 
          self.poptrace.append(min(fP))
        else: 
          self.poptrace.append(max(fP))
      self.result['trace'] = self.trace
      self.result['besttrace'] = self.besttrace
      self.result['poptrace'] = self.poptrace

  def init(self, MAXEVALS, TRACE):
    self.TRACE = TRACE
    self.evals = MAXEVALS
    self.count = 0
    P = self.space.get(self.N)
    fP = self.eval(P)
    self.tracing(P,fP)
    return P, fP

  def run( self, MAXEVALS, TRACE=False ):
    P, fP = self.init(MAXEVALS, TRACE)
    while( not self.stop() ):
      P, fP = self.next(P, fP)
      self.tracing(P,fP)
    if('evals' not in self.result):
      self.result['evals'] = self.count
    return self.result


def experiment(sgoal, problem, MAXEVALS, R=100):
  opt = problem['optimum'] if 'optimum' in problem else None
  r = [sgoal(problem).run(MAXEVALS) for i in range(R)]
  fx = [y['f'] for y in r]
  evals = [y['evals'] for y in r]
  sr = sum([1 if f==opt else 0 for f in fx]) / R
  return fx, evals, sr
