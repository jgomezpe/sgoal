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
from sgoal.util import randbool
from sgoal.util import arity

############### SEARCH SPACE AND  PROBLEM ################

# Simple population generation method 
def simplegetN(N, g):
  return [g() for i in range(N)]

# Search space: a dictionary with
#   'get' : one candidate solution generation method
#   'getN': population generation method. Sets to simplegetN method if no defined
#   'feasible': Method for determining if a candidate solution is feasible. By default every candidate solution is feasible
def SPACE(g, gn = None, feasible = lambda x: True):
  space = {'get':g,'feasible':feasible}
  if(gn==None): gn = lambda N: simplegetN(N,g)
  space['getN'] = gn
  return space

# Problem. Extends the Space dictionary with the following keys
#   'type' : set to 'min' for minimization, set to 'max' for maximization
#   'f' : Wraps the objective function with function eval to record the best progression and traced information
#   'EVALS': Maximum number of function evaluations
#   'minimize': Set to True if the problem is a minimization problem and set to False if a maximization one
#   'pick': Sorts two solutions accoriding to the type of optimization problem (first is the best)
#   'trace': A dictionary with all the information to be traced if required (when set to True in the problem)
#       'f': Evolution (in time) of values of the objective function (calls)
#       'fP':Evolution (in time) of values of the objective function for a population based method
def PROBLEM(type, f, space, EVALS, TRACE=False):
  space['f'] = lambda x: eval(x, f, space)
  space['EVALS'] = EVALS
  space['minimize'] = type=='min'
  if(space['minimize']): space['pick'] = min_pick
  else: space['pick'] = max_pick
  if(TRACE): space['trace'] = {'f':[], 'fP':[], 'best':[]}
  else: space['trace'] = None
  return space

############### GENERIC LIST VARIATION OPERATIONS ################
# Variation 1 to 1: applies a variation and computes the function on the produced candidate solution
def variation11( x, fx, oper, sgoal):  
  y = oper(x)
  fy = sgoal['f'](y)
  return y, fy

# Variation 2 to 1: applies a two parents variation, builds a candidate solution and computes the function on it
def variation21( x, y, fx, fy, oper, sgoal):  
  y1, y2 = oper(x, y)
  if(randbool()): y1, y2 = y2, y1
  fy1 = sgoal['f'](y1)
  return y1, fy1

# Variation 2 to 2: applies a two parents variation, builds two candidate solutions and computes the function on them
def variation22( x, fx, y, fy, oper, sgoal):  
  y1, y2 = oper(x, y)
  if(randbool()): y1, y2 = y2, y1
  fy1 = sgoal['f'](y1)
  fy2 = sgoal['f'](y2)
  return y1, fy1, y2, fy2

# Simple crossover.
def simplexover( x1, x2 ):
  n = len(x1)
  p = rand.randint(1,n-1)
  y1 = x1[0:p] + x2[p:n]
  y2 = x2[0:p] + x1[p:n]
  return y1, y2

# Crossover with mutation.
def xovermutation( x1, x2, xover, mutation, xr=1.0 ):
  if( randbool(xr) ): y1, y2 = xover(x1, x2)
  else: y1, y2 = x1, x2
  if(randbool()): y1, y2 = y2, y1
  return mutation(y1), mutation(y2)

# Transposition
def transposition( x ):
  y = x.copy()
  start = rand.randint(0,len(x)-1)
  end = rand.randint(0,len(x)-1)
  if start>end: start, end = end, start
  while start<end:
    y[start], y[end] = y[end], y[start]
    start += 1
    end -= 1
  return y

############### UTILITY FUNCTIONS ##############
# Returns two candidate solutions according to their function value and minimization problem
# The first solution returned is the best one. 
def min_pick(x, fx, y, fy):
  if(fy <= fx):
    return y, fy, x, fx
  return x, fx, y, fy

# Returns two candidate solutions according to their function value and maximization problem
# The first solution returned is the best one. 
def max_pick(x, fx, y, fy):
  if(fy >= fx):
    return y, fy, x, fx
  return x, fx, y, fy

# Traces a candidate solution and its objective function value
def rec(x, fx, sgoal):
  sgoal['count'] += sgoal['delta']

  if('best' not in sgoal):
    sgoal['best'] ={'x':x, 'f':fx, 'evals':sgoal['count']}
    best = sgoal['best']
  else:
    best = sgoal['best']
    best['x'], best['f'], b, fb = sgoal['pick'](best['x'], best['f'], x, fx)
    if(best['f'] != fb and fx != fb): best['evals'] = sgoal['count']

  if(sgoal['trace'] != None): 
    sgoal['trace']['f'].append(fx)
    sgoal['trace']['best'].append(best['f'])

# Evaluates the objective function in a candidate solution (traces information)
def eval(x, f, sgoal):
  fx = f(x)
  rec(x, fx, sgoal)
  return fx

# Determines if the SGoal can eval the objective function n times
def caneval(sgoal):
  return sgoal['count'] < sgoal['EVALS']

# Determines if the SGoal found the optimum value (if such information is available - for testing purposes)
def optimumfound(sgoal):
  return sgoal['best']['f']==sgoal['optimum']

# Stops the SGoal if cannot eval the objective function anymore or the optimal
# value is found (if available and for testing purposes)
def basicstop(sgoal):
  return not caneval(sgoal) or optimumfound(sgoal)

############### STOCHASTIC GLOBAL OPTIMIZATION ALGORITHM ################
# Extends the Problem dictionary with the following keys
#   'count': Number of function evaluations up to now (when consulted)
#   'delta': Part (or full) function evaluations carried by the last function called (by default set to 1)
#   'optimum': Optimum value if provided. Set to None otherwise
#   'stop': Stopping criteria (predicate without reveiving arguments - usually a lambda of a function receiving the sgoal).
#           By default set to lambda of basicstop
#   'best': Current best solution. A dictionary with the following keys
#       'x': Starting candidate solution
#       'f': value of the objective function on x
#       'evals': Amount of function calls to obtain it
#   'init': The init candidate solution/population method (must be provided by a concrete SGoal)
#   'next': The next candidate solution/population method (must be provided by a concrete SGoal)
#   'start': Optional starting candidate solution. A dictionary with the following keys
#       'x': Starting candidate solution
#       'f': value of the objective function on x
#       'evals': Amount of function calls to obtain it
def SGOAL(problem):
  sgoal = problem
  sgoal['count'] = 0
  sgoal['delta'] = 1
  if('optimum' not in sgoal): sgoal['optimum'] = None
  if('stop' not in sgoal): sgoal['stop'] = lambda : basicstop(sgoal)
  return sgoal

# Runs the SGoal a maximum of MAXEVALS and traces information is desired
def run( sgoal ):
  InitPop, Stop, NextPop, N, trace = sgoal['init'], sgoal['stop'], sgoal['next'], sgoal['N'], sgoal['trace']
  P, fP = InitPop()
  if(N>1): tracepop(fP, trace)
  while( not Stop() ):
    P, fP = NextPop(P, fP)
    if(N>1): 
      sgoal['P'] = P
      sgoal['fP'] = fP
      tracepop(fP, trace)
  return sgoal['best']

##################  SINGLE POINT SGOAL ####################
# Inits a candidate solution
def init(sgoal):
  if('start' in sgoal):
    start = sgoal['start']
    sgoal['best'] = start.copy()
    x = start['x']
    fx = start['f']
    sgoal['count'] = start['evals'] if 'evals' in start else 1
  else:
    x = sgoal['get']()
    fx = sgoal['f'](x)
  return x, fx

# Single point SGoal. Sets the initPopulation to init if not provided
def SPSGoal(problem):
  problem['N'] = 1
  if('init' not in problem): problem['init'] = lambda: init(problem)
  return SGOAL(problem)

########### Variation/Replace Single Point SGOAL ###########
# next method combining a variation and replacement strategies
def next(x, fx, sgoal):
  V, R = sgoal['variation'], sgoal['replace']
  y, fy = V(x, fx)
  return R(x, fx, y, fy)

# Simple replace method: picks the best with neutral mutation
def simplereplace(x, fx, y, fy, sgoal):
  x, fx, y, fy = sgoal['pick'](x, fx, y, fy)
  return x, fx

# Variation/Replace Single Point SGOAL. Extends the SGoal with keys:
#   'variation': Variation operator. Wraps a variation operator if required (when it does not take into account SGOAL info)
#   'replace': Replacement method. Set to simplereplace if not provided
def VRSGoal(problem):
  problem['next'] = lambda x, fx: next(x, fx, problem)
  if('replace' not in problem): problem['replace'] = lambda x,fx, y, fy : simplereplace(x, fx, y, fy, problem)
  if('variation' in problem and arity(problem['variation'])==1): 
    v = problem['variation']
    problem['variation'] = lambda x, fx: variation11(x, fx, v, problem)
  return SPSGoal(problem)

##################  POPULATION SGOAL ####################
# Trace population information
def tracepop(fP, trace):
  if(trace!=None):
    trace['fP'].append(fP.copy())

# Evaluates a population of candidate solutions
def evalPop(P, f): return [f(x) for x in P]

# Inits a population of candidate solutions
def initPop(sgoal):
  f, N, getN = sgoal['f'], sgoal['N'], sgoal['getN']
  if('start' in sgoal):
    start = sgoal['start']
    sgoal['best'] = start.copy()
    x = start['x']
    fx = start['f']
    sgoal['count'] = start['evals'] if 'evals' in start else 1
    if('complement' in sgoal):
      xc = sgoal['complement'](x)
      fxc = f(xc)
      P = getN(N-2)
      fP= evalPop(P, f)
      P.append(xc)
      fP.append(fxc)
    else:  
      P = getN(N-1)
      fP= evalPop(P, f)
    P.append(x)
    fP.append(fx)
  else:
    P = getN(N)
    fP= evalPop(P, f)
  sgoal['P'] = P
  sgoal['fP'] = fP
  return P, fP

# Population Based SGOAL. Sets the population size to 128 and uses initPop as initPopulation method if not provided
def PopSGoal(problem):
  if('N' not in problem): problem['N'] = 128
  if('init' not in problem): problem['init'] = lambda: initPop(problem)
  return SGOAL(problem)


############ EXPERIMENT MODE ###########
# Runs an SGoal in experiment mode. Runs R times the SGoal on the given problem and produces 
#   fx: An array of the best objective function value found by each one of the R runs of the SGoal
#   evals: An array of the objective function evaluations required by the SGoal to achieve such value
#   sr: Success rate. Number of times the SGoal found the optimum value (if available).
def experiment(sgoal, problem, R=100):
  r = []
  for k in range(R):
    p = problem(k)
    opt = p['optimum'] if 'optimum' in p else None
    alg = sgoal(p)
    r.append(run(alg))  
  fx = [y['f'] for y in r]
  evals = [y['evals'] for y in r]
  sr = sum([1 if f==opt else 0 for f in fx]) / R
  return fx, evals, sr