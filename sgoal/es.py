# Evolutionary Strategies, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co and eleonguz@unal.edu.co
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
from sgoal.core import basicInitPop
from sgoal.core import basicStop
from sgoal.binary import bitmutationprob

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# G: Checking evaluations for adapting the mutation rate
# a: Scaling factor of the mutation rate
# mutation: A mutation operation with the possibility of setting its rate
# mr: Initial mutation rate
# x: Initial point
# fx: f value at point x (if provided)
def nextPopR1_5(P, fP, sgoal):
  y = sgoal.mutation(P, sgoal.mr)
  fy = sgoal.evalone(y)
  w = P
  P, fP, y, fy = sgoal.pick( P, fP, y, fy )
  sgoal.I += 1
  if( w==y ): sgoal.Gs += 1
  if( sgoal.I==sgoal.G ):
    p = sgoal.Gs/sgoal.G
    if( p > 0.2 ): sgoal.mr /= sgoal.a
    elif( p < 0.2 ): sgoal.mr *= sgoal.a
    sgoal.Gs = 0
    sgoal.I = 0    
  return P, fP

class Rule_1_5th(SGoal):
  def __init__(self, problem, config, initPop=basicInitPop, stop=basicStop):
    SGoal.__init__(self, problem, nextPopR1_5, initPop, stop)
    self.G = config['G']
    self.a = config['a']
    self.mr = config['mr']
    self.mutation = config['mutation']
    self.I = 1
    self.Gs = 1

def BitArrayR1_5(problem):
  D = problem['space'].D
  return Rule_1_5th(problem, {'mr':1/D, 'mutation': bitmutationprob, 'G':D, 'a':0.9}) 


#def es(f, evals, mutation, marriage, mutation_s, marriage_s, _mu, _lambda, plus, _rho, P, fP=None):
#  if( evals >= len(P) and not fP ):
#    fP = evaluate_population(f, P)
#    evals -= len(P)
  
#  while(evals>0):
#    Q = 