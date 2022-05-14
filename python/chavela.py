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
# promote products derived from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from sgoal import normalize
from sgoal import evaluate
from sgoal import evaluate_population
from sgoal import weighted
from sgoal import arity
from sgoal import pick
from sgoal import tournament
from sgoal import MAXIMIZE
import random as rand

############ RATES OPERATIONS ############
# Init Rates - individual level (M operators)
def init_rates(M):
  rate = []
  for i in range(0,M):
    rate.append(rand.random())
    while rate[i] == 0:
      rate[i] = rand.random()
  return normalize(rate)    

#Init Rates - Population level (A population of N individuals, M operators)
def init_rates_population(N, M):
  return [init_rates(M) for i in range(N)]

# Update rates
def update_rates(h, fc, fp, rate):
  delta = rand.random()
  if ((MAXIMIZE and fc>fp) or (not MAXIMIZE and fc<fp)):
    rate[h] *= (1.0+delta)
  else:
    rate[h] *= (1.0-delta)
  return normalize(rate) 

############ Canonical HAEA ##########
def CHAVELA( f, evals, operators, P, fP=None, rates=None ):
  N = len(P)
  if( not fP and evals>=N ): 
    fP = evaluate_population(f, P)
    evals -= N
  if( not rates ):
    rates = init_rates_population(N,len(operators))
  while(evals>=N):
    Q = []
    fQ = []
    for i in range(N):
      parents = [P[i]]
      h = weighted(rates[i])
      a = arity(operators[h]) 
      if a > 1: 
        idxparents = tournament(fP,a-1)
        for k in idxparents:
          parents.append(P[k])     
        candidates = operators[h](*parents)
        c = candidates[rand.randint(0,len(candidates)-1)]
      else: c = operators[h](parents[0])     
      fc = evaluate(f,c)
      rates[i] = update_rates(h, fc, fP[i], rates[i])
      c, p, fc, fp = pick(P[i], c, fP[i], fc)
      Q.append(c)
      fQ.append(fc)
    evals -= N
    P, fP = Q, fQ
  return P, fP, evals, rates