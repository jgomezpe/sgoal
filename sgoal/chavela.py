# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co  and eleonguz@unal.edu.co
# All rights reserved.
# The Cannonical Hybrid Adaptive Evolutionary Algorithm as proposed by
# Jonatan Gómez and Elizabeth León in: 
# "On the class of hybrid adaptive evolutionary algorithms (chavela)"
# Published in Natural Computing / Issue 3/2021, Print ISSN: 1567-7818 
# Electronic ISSN: 1572-9796, DOI
# https://doi.org/10.1007/s11047-021-09843-5# Publication date 26-02-2021
# Publisher Springer Netherlands
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

from sgoal.core import normalize
from sgoal.core import SGoal
from sgoal.core import weighted
from sgoal.core import arity
from sgoal.core import tournament
import random as rand

#Chavela traced information
RATES = [] # Tracing the evolution of rates
FP = [] # Tracing the population's fitness evolution

# Gets the chavela's traced information
def chavela_trace():
  global RATES,FP
  return FP, RATES
  
############ RATES OPERATIONS ############

# Init Rates - individual level (M operators)
def init_rates(M):
  rate = []
  for i in range(0,M):
    rate.append(rand.random())
    while rate[i] == 0:
      rate[i] = rand.random()
  return normalize(rate)



############ Canonical HAEA: Chavela ##########
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# operators: Genetic operators
# P: Initial population
# fP: Fitness value for each individual of the population, if provided.
# rates: initial rates
# problem: Problem to be solved
# config: GGA parameters
class CHAVELA(SGoal):
  def __init__(self, problem, config):
    SGoal.__init__(self, problem)
    self.N = config['N']
    self.operators = config['operators']
    self.rates = []
    self.result['rates'] = []


  #Init Rates - Population level (A population of N individuals, M operators)
  def init(self, MAXEVALS, TRACE):
    P, fP = SGoal.init(self, MAXEVALS, TRACE)
    M = len(self.operators)
    self.rates = [init_rates(M) for i in range(self.N)]
    return P, fP

  def tracing(self, P, fP):
    SGoal.tracing(self, P, fP)
    if(self.TRACE):
      self.result['rates'].append(self.rates)

  # Determines if the fitness of the second individual is better than the fitness of the first one
  def improves(self, fx, fy):
    return ((not self.minimize and fy>fx) or (self.minimize and fy<fx))

  # Update rates
  def update_rates(self, h, fc, fp, rate):
    delta = rand.random()
    if(self.improves(fp,fc)): 
      rate[h] *= (1.0+delta)
    else: 
      rate[h] *= (1.0-delta)
    return normalize(rate)
    
  def next(self, P, fP):
    Q = []
    fQ = []
    rQ = []
    for i in range(self.N):
      if(self.caneval()):
        parents = [P[i]]
        h = weighted(self.rates[i])
        a = arity(self.operators[h]) 
        if a > 1: 
          idxparents = tournament(fP,a-1, self.minimize)
          for k in idxparents:
            parents.append(P[k])     
          candidates = self.operators[h](*parents)
          c = candidates[rand.randint(0,len(candidates)-1)]
        else: 
          c = self.operators[h](parents[0])     
        fc = self.evalone(c)
        rQ.append(self.update_rates(h, fc, fP[i], self.rates[i]))
        c, fc, p, fp = self.pick(P[i], fP[i], c, fc)
        Q.append(c)
        fQ.append(fc)
      else:
        Q.append(P[i])
        fQ.append(fP[i])
        rQ.append(self.rates[i])
    P, fP, self.rates = Q, fQ, rQ
    return P, fP