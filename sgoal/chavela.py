# The Cannonical Hybrid Adaptive Evolutionary Algorithm as proposed by
# Jonatan Gómez and Elizabeth León in: 
# "On the class of hybrid adaptive evolutionary algorithms (chavela)"
# Published in Natural Computing / Issue 3/2021, Print ISSN: 1567-7818 
# Electronic ISSN: 1572-9796, DOI
# https://doi.org/10.1007/s11047-021-09843-5# Publication date 26-02-2021
# Publisher Springer Netherlands
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

from sgoal.core import caneval
from sgoal.core import SPSGoal
from sgoal.core import PopSGoal
from sgoal.util import normalize
from sgoal.select import weighted
from sgoal.util import arity
from sgoal.select import min_tournament
from sgoal.select import max_tournament
from sgoal.core import initPop
from sgoal.core import init as initcore

import random as rand
  
############ RATES OPERATIONS ############

# Init Rates for an individual (M operators)
def initRates(M):
  return [1/M for i in range(M)]

# Determines if the fitness of the second individual is better than the fitness of the first one
# minimization
def min_improves(fx, fy):
  return fy<fx

# maximization
def max_improves(fx, fy):
  return fy>fx

# Update rates
def updateRates(improves, h, rate):
  delta = rand.random()
  if(improves): 
    rate[h] *= (1.0+delta)
  else: 
    rate[h] *= (1.0-delta)
  return normalize(rate)

# Trace rates
def tracerates(rates, trace):
  if( trace == None): return

  if('rates' not in trace):
    trace['rates'] = []
  r = rates[0].copy()
  for i in range(1, len(rates)):
    for k in range(len(r)):
      r[k] += rates[i][k]
  for k in range(len(r)):
    r[k] /= len(rates)
  trace['rates'].append(r)

def apply( oper, x, fx, sgoal ):
  P, fP, selection = sgoal['P'], sgoal['fP'], sgoal['selection']
  a = arity(oper) 
  if a > 1: 
    parents = [x]
    idxparents = selection(fP,a-1)
    for k in idxparents:
      parents.append(P[k])
    c = oper(*parents)
    c = c[rand.randint(0,len(c)-1)]
  else: 
    c = oper(x)
  fc = sgoal['f'](c)
  return c, fc

# CHAVELA Inits Population: Generates a population (N individuals) following an inner 
# population generation process (by default the space population generation)
# and associates operator rates to each individual (M operators)
def init(sgoal):
  N = sgoal['N']
  P, fP = sgoal['innerInit'](sgoal)
  M = len(sgoal['operators'])
  sgoal['rates'] = [initRates(M) for i in range(N)]
  tracerates(sgoal['rates'], sgoal['trace'])
  return P, fP

# CHAVELA next population method
def next(P, fP, sgoal):
  improves, operators, pick, N, rates = sgoal['improves'], sgoal['operators'], sgoal['pick'], sgoal['N'], sgoal['rates']
  Q = []
  fQ = []
  rQ = []
  for i in range(N):
    if(caneval(sgoal)):
      h = weighted(rates[i])
      c, fc = operators[h](P[i], fP[i])
      rQ.append(updateRates(improves(fP[i], fc), h, rates[i]))
      c, fc, p, fp = pick(P[i], fP[i], c, fc)
      Q.append(c)
      fQ.append(fc)
    else:
      Q.append(P[i])
      fQ.append(fP[i])
      rQ.append(rates[i])
  P, fP, sgoal['rates'] = Q, fQ, rQ
  tracerates(sgoal['rates'], sgoal['trace'])
  return P, fP

############ Generic Canonical HAEA: Chavela ##########
# Extends the Population SGOAL with the following keys:
#   'rates': operator's rates one by each candidate solution in the population
#   'operators': A list with the variation operators used by CHAVELA 
#   'innerInit': Inner Init Population method by default set to initPop from sgoal.core
def CHAVELA_T(problem):
  if(problem['minimize']): 
    problem['improves'] = min_improves
    if('selection' not in problem): problem['selection'] = min_tournament
  else: 
    problem['improves'] = max_improves
    problem['selection'] = max_tournament
  if('init' in problem): problem['innerInit'] = problem['init']
  else: problem['innerInit'] = initPop
  problem['init'] = lambda : init(problem)
  problem['next'] = lambda P, fP: next(P, fP, problem)
  problem['rates'] = []
  return PopSGoal(problem)

############### Single Point CHAVELA ############
# generates one candidate solution following an inner single point generation strategy
def init1(sgoal):
  x, fx = sgoal['innerInit'](sgoal)
  M = len(sgoal['operators'])
  sgoal['rates'] = initRates(M)
  return x, fx

# CHAVELA1 next individual method
def next1(x, fx, sgoal):
  f, minimize, operators, pick, rates = sgoal['f'], sgoal['minimize'], sgoal['operators'], sgoal['pick'], sgoal['rates']
  if(minimize): improves = min_improves
  else: improves = max_improves
  if(not caneval(sgoal)): return x, fx
  h = weighted(rates)
  c = operators[h](x)
  fc = f(c)
  sgoal['rates'] = updateRates(improves(fx,fc), h, rates)
  x, fx, c, fc = pick(x, fx, c, fc)
  if('tracerates' in sgoal and sgoal['tracerates']): 
    sgoal['tracerates'](sgoal)    
  return x, fx


############ Single point CHAVELA ##########
# Extends the Single Point SGOAL with the following keys:
#   'rates': operator's rates one for the candidate solution
#   'operators': A list with the variation operators used by CHAVELA 
#   'innerInit': Inner Init candidate solution method by default set to init from sgoal.core
def CHAVELA1_T(problem):
  if('init' in problem): problem['innerInit'] = problem['init']
  else: problem['innerInit'] = initcore
  problem['init'] = lambda: init1(problem)
  problem['next'] = lambda x, fx: next1(x, fx, problem)
  problem['rates'] = []
  return SPSGoal(problem)