# One dimensional Real space definitions
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

########### Interval Space [min,max] ########### 
# Creates a real vector interval/hyperrectangle (according to D and values min and max) as search space
def Interval(min, max, g=None, gn=None, feasible=None):
  length = max - min
  if(g==None): g = lambda: min + rand.random()*length
  if(gn==None): gn = lambda N: simplegetN(N,g)
  if(feasible==None): feasible = lambda x: min <= x <= max
  space = SPACE(g, gn, feasible)
  space['interval'] = [min, max, length]
  space['D'] = 1
  return space

# Gaussian mutation
def intervalGaussian(x, sigma, feasible):
  y = x + rand.gauss(sigma)
  return y if feasible(y) else x

def gaussiansigma(sigma, feasible): return lambda x: intervalGaussian(x, sigma, feasible)

def gaussian(sgoal): 
  feasible = sgoal['feasible']
  sigma = sgoal['interva'][2]/100.0
  return lambda x: intervalGaussian(x, sigma, feasible)

#################### TEST FUNCTIONS ##############
# Sphere function
def sphere_1(x):
   return x*x

# Rastrigin function as proposed by Rastrigin, L. A. in "Systems of extremal control." Mir, Moscow (1974).
def rastrigin_1( x ):
	return x*x - 10.0*math.cos(2.0*math.pi*x)

# Schwefel Function
def schwefel_1( x ):
	return -x * math.sin(math.sqrt(abs(x)))

# Griewangk function
def griewank_1( x ):
  sum = x*x/4000.0
  prod = math.cos(x)
  return (1.0 + sum - prod)

def ackley_1( x ):
  e = math.exp(1)
  a = 20
  b = 0.2
  c = 2*math.pi
  s = a + e
  p = x*x
  q = math.cos(c*x)
  return s-a*math.exp(-b*p**0.5) - math.exp(q)

##################### TEST PROBLEMS ####################
def RealTestProblem(f, D, EVALS, TRACE=False):
  if(f=='Rastrigin'): problem = PROBLEM('min', rastrigin_1, Interval(-5.12, 5.12), EVALS, TRACE)
  elif(f=='Schwefel'): problem = PROBLEM('min', schwefel_1, Interval(-500.0, 500.0), EVALS, TRACE)
  elif(f=='Griewank'): problem = PROBLEM('min', griewank_1, Interval(-600.0, 600.0), EVALS, TRACE)
  elif(f=='Ackley'): problem = PROBLEM('min', ackley_1, Interval(-32.768, 32.768), EVALS, TRACE)
  elif(f=='Sphere'): problem = PROBLEM('min', sphere_1, Interval(-5.12, 5.12), EVALS, TRACE)
  else: problem = PROBLEM('min', sphere_1, Interval(-5.12, 5.12), EVALS, TRACE)
  problem['optimum'] = 0.0
  return problem
