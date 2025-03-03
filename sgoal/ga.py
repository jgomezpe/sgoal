# Genetic Algorithms as proposed by 
# J. H. Holland, Adaptation in Natural and Artificial Systems. 
# The University of Michigan Press, 1975.
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
from sgoal.core import randbool
from sgoal.core import PopSGoal
from sgoal.core import variation11
from sgoal.core import variation22
from sgoal.core import simplexover
from sgoal.core import xovermutation
from sgoal.select import min_tournament
from sgoal.select import max_tournament
from sgoal.util import arity
import random as rand


############### Genetic Algorithm - GA ################
# problem: Problem to solve
# config: GA Configuration
#   selection: Selection Mechanism
#   mutation: Mutation operator
#   xr: Crossover rate
#   xover: Crossover operator
#   N: Population's size (we set to the closest higher even number)
# nextPop: Next population generation method
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
def GA(problem):
  if('N' in problem and problem['N'] % 2 != 0): problem['N'] += 1
  if( 'selection' not in problem ): 
    if(problem['minimize']): problem['selection'] = min_tournament
    else: problem['selection'] = max_tournament 
  if('nextpair' not in problem):
    if( 'xr' not in problem ):  problem['xr'] = 0.7
    mut = problem['mutation']
    if( 'xover' not in problem ): xover = simplexover
    else: xover = problem['xover']
    problem['nextpair'] = lambda x, y : xovermutation(x, y, xover, mut, problem['xr'])
  np = problem['nextpair']
  if(arity(np)==2):
    problem['nextpair'] = lambda x, fx, y, fy: variation22(x, fx, y, fy, np, problem)
  elif(arity(np)==5):
    problem['nextpair'] = lambda x, fx, y, fy: np(x, fx, y, fy, problem)
    
  return PopSGoal(problem)

############### Generational Genetic Algorithm - GGA ################
def nextGGA(P, fP, sgoal):
  N, selection, nextpair = sgoal['N'], sgoal['selection'], sgoal['nextpair']
  idx1, idx2 = selection(fP, 2)
  Q = []
  fQ = []
  for i in range(N//2):
    if(caneval(sgoal)):
      a, fa, b, fb = nextpair(P[idx1], fP[idx1], P[idx2], fP[idx2])
      Q.append(a)
      fQ.append(fa)
      Q.append(b)
      fQ.append(fb)
  return Q, fQ

############### Generic Generational Genetic Algorithm - GGA ################
# problem: Problem to solve
# config: GA Configuration
#   selection: Selection Mechanism
#   mutation: Mutation operator
#   xr: Crossover rate
#   xover: Crossover operator
#   N: Population's size (we set to the closest higher even number)
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
def GGA_T(problem):
  problem['next'] = lambda P, fP: nextGGA(P, fP, problem)
  return GA(problem)

############### Steady State Genetic Algorithm - GGA ################
def nextSSGA(P, fP, sgoal):
  N, nextpair, pick = sgoal['N'], sgoal['nextpair'], sgoal['pick']
  for i in range(N//2):
    if(caneval(sgoal)):
      idx1, idx2 = rand.randint(0, N-1), rand.randint(0, N-1)
      a, fa, b, fb = nextpair(P[idx1], fP[idx1], P[idx2], fP[idx2])
      P[idx1], fP[idx1], a, fa = pick(P[idx1], fP[idx1], a, fa)
      P[idx2], fP[idx2], b, fb = pick(P[idx2], fP[idx2], b, fb)
  return P, fP

############### Generic Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
# config: GA Configuration
#   selection: Selection Mechanism
#   mutation: Mutation operator
#   xr: Crossover rate
#   xover: Crossover operator
#   N: Population's size (we set to the closest higher even number)
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
def SSGA_T(problem):
  problem['next'] = lambda P, fP: nextSSGA(P, fP, problem)
  return GA(problem)
