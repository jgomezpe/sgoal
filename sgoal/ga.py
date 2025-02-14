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
from sgoal.core import SGoal
from sgoal.core import randbool
from sgoal.core import best
from sgoal.core import basicInitPop
from sgoal.core import basicStop
from sgoal.core import tournament
from sgoal.core import simplexover
from sgoal.binary import bitmutation
from sgoal.real import lambdaGaussianMutation
from sgoal.gabo import initPopGABO2N

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
class GA(SGoal):
  def __init__(self, problem, config, nextPop, initPop=basicInitPop, stop=basicStop):
    SGoal.__init__(self, problem, nextPop, initPop, stop)
    self.selection = config['selection']
    self.mutation = config['mutation']
    self.xr = config['xr']
    self.xover = config['xover']
    self.N = config['N']
    if(self.N%2==1): 
      self.N+=1
    self.poptrace = []


############### Generational Genetic Algorithm - GGA ################
def nextPopGGA(P, fP, sgoal):
  Q = []
  fQ = []
  for i in range(0,sgoal.N,2):
    idx1, idx2 = sgoal.selection(fP, 2, sgoal.minimize)
    if sgoal.caneval() and randbool(sgoal.xr):
      a, b = sgoal.xover(P[idx1], P[idx2])
      a = sgoal.mutation(a)
      b = sgoal.mutation(b)
      Q.append(a)
      fQ.append(sgoal.evalone(a))
      if(sgoal.caneval()):
        Q.append(b)
        fQ.append(sgoal.evalone(b) )
    else:
      Q.append( P[idx1] )
      Q.append( P[idx2] )
      fQ.append(fP[idx1] )
      fQ.append(fP[idx2] )
  P = Q
  fP = fQ
  return P, fP

############### Generational Genetic Algorithm - GGA ################
# problem: Problem to solve
# config: GA Configuration
#   selection: Selection Mechanism
#   mutation: Mutation operator
#   xr: Crossover rate
#   xover: Crossover operator
#   N: Population's size (we set to the closest higher even number)
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
def GGA(problem, config, initPop=basicInitPop, stop=basicStop):
  return GA(problem, config, nextPopGGA, initPop, stop)

############### Steady State Genetic Algorithm - GGA ################
def nextPopSSGA(P, fP, sgoal):
  i=0
  while(i<sgoal.N and sgoal.caneval()):
    idx1, idx2 = sgoal.selection(fP, 2, sgoal.minimize)
    p1, p2 = P[idx1], P[idx2]
    a, b = sgoal.xover(p1, p2)
    a, b = sgoal.mutation(a), sgoal.mutation(b)
    fa = sgoal.evalone(a)
    k = best(fP, not sgoal.minimize)
    P[k] = a
    fP[k] = fa
    if(sgoal.caneval()):
      k = best(fP, not sgoal.minimize)
      fb= sgoal.evalone(b)
      P[k] = b
      fP[k] = fb
    i += 2
  return P, fP

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
# config: GA Configuration
#   selection: Selection Mechanism
#   mutation: Mutation operator
#   xr: Crossover rate
#   xover: Crossover operator
#   N: Population's size (we set to the closest higher even number)
# initPop: Process for generating the initial population (by default uses the BitArraySpace generation method)
# stop: Stopping criteria (by default uses the basic stopping criteria)
def SSGA(problem, config, initPop=basicInitPop, stop=basicStop):
  return GA(problem, config, nextPopSSGA, initPop, stop)

############# Binary versions ##############
# Basic Genetic Algorithm configuration
#   selection: tournament selection
#   xover: simple crossover
#   xr; Crossover rate set to 0.7
#   mutation: bitmutation
#   N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the bitarray length
def BitArrayGAconfig(D, N=-1):
  if(N<0): N = max(D//2,100)
  return {'selection': tournament, 'xover': simplexover, 'xr':0.7, 'mutation':bitmutation, 'N':N}
  
############### Generational Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def BitArrayGGA(problem, N=-1):
  return GGA(problem, BitArrayGAconfig(problem['space'].D, N))

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def BitArraySSGA(problem, N=-1):
  return SSGA(problem, BitArrayGAconfig(problem['space'].D, N))

############### Generational Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def GGA_G(problem, N=-1):
  return GGA(problem, BitArrayGAconfig(problem['space'].D, N), initPopGABO2N)

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def SSGA_G(problem, N=-1):
  return SSGA(problem, BitArrayGAconfig(problem['space'].D, N), initPopGABO2N)

############# Real versions ##############
# Basic Genetic Algorithm configuration for real encoding
#   selection: tournament selection
#   xover: simple crossover
#   xr; Crossover rate set to 0.7
#   mutation: gaussian mutation with sigma=0.2
#   N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def RealGAconfig(space, N=-1):
  D = len(space.min) if isinstance(space.min, list) else 1
  if(N<0): N = max(D//2,100)
  return {'selection': tournament, 'xover': simplexover, 'xr':0.7, 'mutation':lambdaGaussianMutation(0.2, space), 'N':N}
  
############### Generational Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def RealGGA(problem, N=-1):
  return GGA(problem, RealGAconfig(problem['space'], N))

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
# N: Population's size. If N<0 then the population's size is set to max{D//2, 100}, with D the space dimension
def RealSSGA(problem, N=-1):
  return SSGA(problem, RealGAconfig(problem['space'], N))