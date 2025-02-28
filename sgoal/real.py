# Real vector search space definitions
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

import math
import random as rand
from sgoal.core import randbool
from sgoal.core import SPACE
from sgoal.core import PROBLEM
from sgoal.core import simplegetN
from sgoal.core import VRSGoal
from sgoal.core import variation
from sgoal.core import transposition
from sgoal.core import simplexover
from sgoal.real1 import rastrigin_1
from sgoal.real1 import schwefel_1
from sgoal.es import Rule1_5_T
from sgoal.ga import SSGA_T
from sgoal.ga import GGA_T
from sgoal.chavela import CHAVELA_T


########### MyperRectangle/HyperCube Space [min,max] with min and max n-dimensional real vectors ###########
# Checks if a real vector x is in the hyperrectangle [min,max]
def feasibleRn(min, max, x):
  for i in range(len(x)):
    if(not (min[i] <= x[i] <= max[i])):
        return False
  return True

# Creates a real vector interval/hyperrectangle (according to values min and max) as search space
def HyperRectangle(min, max, g=None, gn=None, feasible=None):
  length = [max[i]-min[i] for i in range(len(min))]
  if(g==None): g = lambda: [min[i] + rand.random()*length[i] for i in range(len(min))]
  if(gn==None): gn = lambda N: simplegetN(N,g)
  if(feasible==None): feasible = lambda x: feasibleRn(min, max, x)
  space = SPACE(g, gn, feasible)
  space['D'] = len(min)
  space['hyperrectangle'] = [min, max, length]
  return space

def HyperCube(min, max, D, g=None, gn=None, feasible=None):
  return HyperRectangle([min for i in range(D)], [max for i in range(D)], g, gn, feasible)

############# Variations ##############
# N-Dimensional Gaussian mutation
def hyperGaussianSigmaProb(x, sigma, feasible, p):
  y = x.copy()
  for i in range(len(y)):
     if(randbool(p)):
        y[i] += rand.gauss(sigma[i])
  return y if feasible(y) else x

def hypergaussiansigma(x, sigma, feasible):
  return hyperGaussianSigmaProb(x, sigma, feasible, 1.0/len(x))

def hypergaussian(sigma, feasible): return lambda x: hypergaussiansigma(x, sigma, feasible)

def hypergaussianmutation(problem):
  if( 'hyperrectangle' in problem):
    L = problem['hyperrectangle'][2].copy()
    sigma = [s/100 for s in L]
  else:
    D = problem['D']
    sigma = [0.2 for i in range(D)]
  problem['sigma'] = sigma
  feasible = problem['feasible']
  return hypergaussian(sigma, feasible)

def gaussianmutation(sgoal):
  if('mutation' not in sgoal): sgoal['mutation'] = hypergaussianmutation(sgoal)
  return sgoal

##################### SGOALs ###########################
# Classical Hill Climbing Algorithm for Real problems. Uses Gaussian mutation with sigma=0.2 as variation operator
# problem: Problem to solve
def HC(problem): 
  if( 'variation' not in problem ): problem['variation'] = hypergaussianmutation(problem)
  return VRSGoal(problem) 

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
# Scale the current sigma parameter
def scalesigma( problem ):
  v = problem['parameter']
  sigma = problem['sigma']
  sigma = [s*v for s in sigma]
  problem['sigma'] = sigma
  problem['variation'] = lambda x, fx: variation(x, fx, hypergaussianmutation(problem), problem)

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule
def Rule1_5(problem):
  D = problem['D']
  if( 'parameter' not in problem ): problem['parameter'] = 1
  if( 'variation' not in problem ): 
    problem['scaleparameter'] = lambda : scalesigma(problem)
    problem['variation'] = hypergaussianmutation(problem)
  if( 'G' not in problem ): problem['G'] = D
  return Rule1_5_T(problem)

############## Generational Genetic Algorithm - SSGA ################
def GGA(problem):
  return GGA_T(gaussianmutation(problem))

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
def SSGA(problem):
  return SSGA_T(gaussianmutation(problem))

# Standard CHAVELA for Real problems. Uses gaussianmutation, simplexover, and transposition as operators
def CHAVELA(problem):
  if( 'operators' not in problem ): problem['operators'] = [hypergaussianmutation(problem), simplexover, transposition]
  return CHAVELA_T(problem) 

#################### TEST FUNCTIONS ##############
# Sphere function
def sphere(x):
   s = 0.0
   for y in x:
      s += y*y
   return s

# Rastrigin function as proposed by Rastrigin, L. A. in "Systems of extremal control." Mir, Moscow (1974).
def rastrigin( x ):
  f = 0.0
  for c in x:
    f += rastrigin_1(c)
  return 10.0*len(x) + f

# Schwefel Function
def schwefel( x ):
  f = 0.0
  for c in x:
    f += schwefel_1(c)
  return (418.9829101*len(x) + f)

# Griewangk function
def griewank( x ):
  sum = 0.0
  prod = 1.0
  for i in range(len(x)):
    sum += x[i]*x[i]/4000.0
    prod *= math.cos(x[i]/math.sqrt(i+1.0))
  return (1.0 + sum - prod)

# Rosenbrock Saddle Function
def rosenbrock_saddle_2( x1, x2 ):
	y = x1*x1 - x2
	return (100.0*y*y + (1.0-x1)*(1.0-x1))

def rosenbrock_saddle( x ):
  f = 0.0
  for i in range(len(x)-1):
    f += rosenbrock_saddle_2( x[i], x[i+1] )
  return f

def ackley( x ):
  D = len(x)
  e = math.exp(1)
  a = 20
  b = 0.2
  c = 2*math.pi
  s = a + e
  p = 0
  q = 0
  for y in x:
    p += y*y
    q += math.cos(c*y)
  return s-a*math.exp(-b*(p/D)**0.5) - math.exp(q/D)

##################### TEST PROBLEMS ####################
def TestProblem(f, D, EVALS, TRACE=False):
  if(f=='Rastrigin'): problem = PROBLEM('min', rastrigin, HyperCube(-5.12, 5.12, D), EVALS, TRACE)
  elif(f=='Schwefel'): problem = PROBLEM('min', schwefel, HyperCube(-500.0, 500.0, D), EVALS, TRACE)
  elif(f=='Griewank'): problem = PROBLEM('min', griewank, HyperCube(-600.0, 600.0, D), EVALS, TRACE)
  elif(f=='Rosenbrock'): problem = PROBLEM('min', rosenbrock_saddle, HyperCube(-2.048, 2.048, D), EVALS, TRACE)
  elif(f=='Ackley'): problem = PROBLEM('min', ackley, HyperCube(-32.768, 32.768, D), EVALS, TRACE)
  elif(f=='Sphere'): problem = PROBLEM('min', sphere, HyperCube(-5.12, 5.12, D), EVALS, TRACE)
  else: problem = PROBLEM('min', sphere, HyperCube(-5.12, 5.12, D), EVALS, TRACE)
  problem['optimum'] = 0.0
  return problem