# Stochastic Global Optimization Algorthms (SGoals) library core functions.
# A formal description is provided by Jonatan Gomez in
# "Stochastic global optimization algorithms: A systematic formal approach",
# Information Sciences, Volume 472, 2019, Pages 53-76, ISSN 0020-0255,
# https://doi.org/10.1016/j.ins.2018.09.021.
# (https://www.sciencedirect.com/science/article/pii/S0020025517305248)
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

# Tournament selection. Selects 1 individual.
# Picks 4 individuals at random and returns the best one
def tournament1(quality, MINIMIZE=True):
  m = 4 #Tournament's size
  candidate = uniform(quality, m, MINIMIZE)
  x = 0
  for k in range(1,m):
    if (MINIMIZE and quality[candidate[x]] >= quality[candidate[k]]) or (not MINIMIZE and quality[candidate[x]] <= quality[candidate[k]]):
      x = k
  return candidate[x]
  
# Tournament selection. Selects N individuals.
# Uses tournament1 for each individual to be selected
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

# Adjust function values (quality arrays) to real quality measures (q_i > 0)
def adjustquality(quality, MINIMIZE=True):
  if(MINIMIZE):
    quality = [-q for q in quality]
  else:
    quality = quality.copy()
  m = min(quality)
  i=0
  n = len(quality)
  while(i<n and m==quality[i]):
    i+=1
  if(i==n):
    return [1 for i in range(n)]
  m2 = quality[i]
  i+=1
  while(i<n):
    if(m < quality[i] and m2>quality[i]):
      m2 = quality[i]
    i+=1
  d = m2 - m
  return [q - m + d for q in quality]

# Roulette wheel selection. Selects N individuals.
def roulette(quality, N, MINIMIZE=True):
  weight = normalize(adjustquality(quality, MINIMIZE))
  p = normalize(weight)
  return [weighted(p) for i in range(N)]

############### SEARCH SPACE ################
class Space:
  # Gets one point in the search space
  def getone(self):
    return []

  # Determines if a candidate solution is feasible
  def feasible(self, x):
    return True

  # Gets N points in the search space
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
# Stops the SGoal if cannot eval the objective function anymore or the optimal
# value is found (if available and for testing purposes)
def basicStop(sgoal):
  return not sgoal.caneval() or sgoal.optimumfound()

# Basic population initialization
def basicInitPop(sgoal):
  P = sgoal.space.get(sgoal.N)
  fP = sgoal.eval(P)
  sgoal.tracing(P,fP)
  return P, fP

#Stochastic Global Optimization Algorithm
# problem: Problem to solve
#   f: Objective function
#   space: Solution space 
#   minimize: If minimizing (True) or maximizing (False)
# nextPop: Next population generation method
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
class SGoal:
  def __init__(self, problem, nextPop, initPop=basicInitPop, stop=basicStop):
    self.f = problem['f']
    self.space = problem['space']
    self.minimize = problem['type'].lower()=='min'
    self.optimum = problem['optimum'] if 'optimum' in problem else None
    self.count = 0
    self.result = {}
    self.trace = []
    self.poptrace = []
    self.besttrace = []
    self.N = 1
    self.initPop = initPop
    self.nextPop = nextPop
    self.stop = stop
    self.evals = 1
    self.TRACE = False

  # Runs the SGoal a maximum of MAXEVALS and tracing information is desired
  def run( self, MAXEVALS, TRACE=False ):
    self.evals = MAXEVALS
    self.TRACE = TRACE
    P, fP = self.initPop(self)
    while( not self.stop(self) ):
      P, fP = self.nextPop(P, fP, self)
      if(TRACE):
        self.tracing(P,fP)
    if('evals' not in self.result):
      self.result['evals'] = self.count
    return self.result

  # Determines if the SGoal can eval the objective function n times
  def caneval(self, n=1):
    return self.count+n <= self.evals
  
  # Determines if the SGoal found the optimum value (if such information is available - for testing purposes)
  def optimumfound(self):
    return self.result['f']==self.optimum
  
  # Trace method. Traces best value found, population information, and computed objective function values
  def tracing(self, P, fP):  
    if(self.N > 1):
      if(self.minimize): 
        self.poptrace.append(min(fP))
      else: 
        self.poptrace.append(max(fP))
    self.result['trace'] = self.trace
    self.result['besttrace'] = self.besttrace
    self.result['poptrace'] = self.poptrace

  # Returns two candidate solutions according to their function value and type of optimization problem
  # The first solution returned is the best one. 
  def pick(self, x, fx, y, fy):
    if((self.minimize and fy <= fx) or (not self.minimize and fy >= fx)):
        return y, fy, x, fx
    return x, fx, y, fy

  # Evaluates the objective function in a candidate solution (traces information)
  def evalone(self, x):
    self.count+=1
    fx = self.f(x)
    if('x' not in self.result):
      self.result['x'] = x
      self.result['f'] = fx
    else:
      self.result['x'], self.result['f'], b, fb = self.pick(x, fx, self.result['x'], self.result['f'])
    if(fx==self.optimum):
      self.result['evals'] = self.count
    if(self.TRACE):
      self.trace.append(fx)
      self.besttrace.append(self.result['f'])
    return fx

  # Evaluates the objective function in a population (traces information accordingly)
  def evalpop(self, P):
    i=0
    fP = []
    while(self.caneval() and i<len(P)):
      fP.append(self.evalone(P[i]))
      i+=1
    return fP

  # Evals a candidate solution or a population.
  def eval(self, x):
    if(self.N==1):
      return self.evalone(x)
    return self.evalpop(x)

# Runs an SGoal in experiment mode. Runs R times the SGoal on the given problem and produces 
# fx: An array of the best objective function value found by each one of the R runs of the SGoal
# evals: An array of the objective function evaluations required by the SGoal to achieve such value
# sr: Success rate. Number of times the SGoal found the optimum value (if available).
def experiment(sgoal, problem, MAXEVALS, R=100):
  opt = problem['optimum'] if 'optimum' in problem else None
  r = [sgoal(problem).run(MAXEVALS) for i in range(R)]
  fx = [y['f'] for y in r]
  evals = [y['evals'] for y in r]
  sr = sum([1 if f==opt else 0 for f in fx]) / R
  return fx, evals, sr