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

from sgoal.core import SGoal
from sgoal.core import normalize
from sgoal.core import transposition
from sgoal.core import weighted
from sgoal.core import arity
from sgoal.core import tournament
from sgoal.core import basicInitPop
from sgoal.core import basicStop
from sgoal.core import tournament
from sgoal.core import simplexover
from sgoal.binary import bitmutation

import random as rand
  
############ RATES OPERATIONS ############

# Init Rates for an individual (M operators)
def initRates(M):
  rate = []
  for i in range(0,M):
    rate.append(rand.random())
    while rate[i] == 0:
      rate[i] = rand.random()
  return normalize(rate)

# CHAVELA Inits Population: Generates a population (N individuals) following an inner 
# population generation process (by default the space population generation)
# and associates operator rates to each individual (M operators)
def initPopCHAVELA(sgoal):
  P, fP = sgoal.innerInit(sgoal)
  M = len(sgoal.operators)
  sgoal.rates = [initRates(M) for i in range(sgoal.N)]
  return P, fP

# CHAVELA next population method
def nextPopCHAVELA(P, fP, sgoal):
  Q = []
  fQ = []
  rQ = []
  for i in range(sgoal.N):
    if(sgoal.caneval()):
      parents = [P[i]]
      h = weighted(sgoal.rates[i])
      a = arity(sgoal.operators[h]) 
      if a > 1: 
        idxparents = tournament(fP,a-1, sgoal.minimize)
        for k in idxparents:
          parents.append(P[k])     
        candidates = sgoal.operators[h](*parents)
        c = candidates[rand.randint(0,len(candidates)-1)]
      else: 
        c = sgoal.operators[h](parents[0])     
      fc = sgoal.evalone(c)
      rQ.append(sgoal.updateRates(h, fc, fP[i], sgoal.rates[i]))
      c, fc, p, fp = sgoal.pick(P[i], fP[i], c, fc)
      Q.append(c)
      fQ.append(fc)
    else:
      Q.append(P[i])
      fQ.append(fP[i])
      rQ.append(sgoal.rates[i])
  P, fP, sgoal.rates = Q, fQ, rQ
  return P, fP

############ Canonical HAEA: Chavela ##########
# problem: Problem to solve
# config: CHAVELA parameters
#   N: Number of individuals
#   operators: A list with the variation operators used by CHAVELA 
# initPop: Process for generting the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
class CHAVELA(SGoal):
  def __init__(self, problem, config, initPop=basicInitPop, stop=basicStop):
    SGoal.__init__(self, problem, nextPopCHAVELA, initPopCHAVELA, stop)
    self.N = config['N']
    self.operators = config['operators']
    self.rates = []
    self.result['rates'] = []
    self.innerInit = initPop

  # Trace method. Adds rates evolution to the tracing object
  def tracing(self, P, fP):
    SGoal.tracing(self, P, fP)
    if(self.TRACE):
      self.result['rates'].append(self.rates)

  # Determines if the fitness of the second individual is better than the fitness of the first one
  def improves(self, fx, fy):
    return ((not self.minimize and fy>fx) or (self.minimize and fy<fx))

  # Update rates
  def updateRates(self, h, fc, fp, rate):
    delta = rand.random()
    if(self.improves(fp,fc)): 
      rate[h] *= (1.0+delta)
    else: 
      rate[h] *= (1.0-delta)
    return normalize(rate)
    

def BitArrayCHAVELA(problem):
  return CHAVELA(problem, {'operators':[bitmutation, simplexover, transposition], 'N':problem['space'].D//2}) 
