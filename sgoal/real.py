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
from sgoal.binary import Binary

# Checks if a real value x is in the interval [min,max]
def feasibleR(min, max, x):
  return min <= x <= max

# Checks if a real vector x is in the hyperrectangle [min,max]
def feasibleRn(min, max, x):
  for i in range(len(x)):
    if(not (min[i] <= x[i] <= max[i])):
        return False
  return True

# Gets a real in the interval [min,max]
def getR(min, length):
  return min + rand.random()*length

# Gets a real vector in the hyperrectangle [min,max]
def getRn(min, length):
  return [min[i] + rand.random()*length[i] for i in range(len(min))]

# Creates a real vector interval/hyperrectangle (according to D and values min and max) as search space
def RealSpace(min, max, D=0):
  if(D>1):
    min = [min for i in range(D)]
    max = [max for i in range(D)]
  if(isinstance(min , list)):
    length = [max[i]-min[i] for i in range(len(min))]
    g = lambda: getRn(min, length)
    space = SPACE(g, lambda N: simplegetN(N,g), lambda x: feasibleRn(min, max, x))
    space['D'] = len(min)
  else:
    length = max - min
    g = lambda: getR(min, length)
    space = SPACE(g, lambda N: simplegetN(N,g), lambda x: feasibleR(min, max, x))
    space['D'] = 1
  space['RSMIN'] = min
  space['RSLENGTH'] = length
  return space

# Gaussian number generation
def gaussian(sigma=1.0): return rand.gauss(sigma)

# Power law number generation
def powerlaw(a=-2.0):
  if a==-2.0: return 1/(1-rand.random())
  else: return (1-rand.random())**(1/(a+1))

############### VARIATION OPERATIONS ################
# Applies a mutation to each real component with probability p
def intensity_mutation(x, mutation, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): y[i] = mutation(y[i]) 
  return y

# Gaussian Mutation for real numbers
def gaussianMutationR( x, sigma, feasible ):
  y = x + gaussian(sigma)
  return y if feasible(y) else x

# Gaussian Mutation for real vectors
def gaussianMutationRn( x, sigma, feasible ):
  y = x.copy()
  i = rand.randint(0,len(x)-1)
  y[i] += gaussian(sigma) # = [z + gaussian(sigma) for z in x]
  return y if feasible(y) else x

# Gaussian Mutation for real vectors spaces
def gaussianMutation(sgoal):
  D = sgoal['D']
  if( D==1 ): return lambda x: gaussianMutationR(x, 0.2, sgoal['feasible'])
  return lambda x: gaussianMutationRn(x, 0.2, sgoal['feasible'])


##################### REAL TO BINARY #####################
# Grows a binary representation to real vector representation
def grow1(x, min, length, SIZE, i=0):
  s = 0
  p = 1
  start = i*SIZE
  for k in range(start, start+SIZE):
    if( x[k] == 1):
      s += p
    p *= 2
  return min + (s/p)*length

def grow(x, min, length, SIZE=16):
  return [grow1(x, min[i], length[i], SIZE, i) for i in range(len(x)//SIZE)]

def Real2Binary(min, max, D=0, BITSIZE=32):
  if(D>1):
    min = [min for i in range(D)]
    max = [max for i in range(D)]
    D *= BITSIZE
  else:
    D = BITSIZE
  space = Binary(D)
  if(isinstance(min , list)):
    length = [max[i]-min[i] for i in range(len(min))]
    space['feasible'] = lambda x: feasibleRn(min, max, x)
    space['grow'] = lambda x: grow(x, min, length, BITSIZE)
  else:
    length = max - min
    space['feasible'] = lambda x: feasibleR(min, max, x)
    space['grow'] = lambda x: grow1(x, min, length, BITSIZE)
  return space

def Real2BinaryPROBLEM(type, f, space, EVALS, TRACE=False):
  return PROBLEM(type, lambda x: f(space['grow'](x)), space, EVALS, TRACE)

#################### TEST FUNCTIONS ##############
# Sphere function
def sphere_1(x):
   return x*x

def sphere(x):
   s = 0.0
   for y in x:
      s += y*y
   return s

# Rastrigin function as proposed by Rastrigin, L. A. in "Systems of extremal control." Mir, Moscow (1974).
def rastrigin_1( x ):
	return x*x - 10.0*math.cos(2.0*math.pi*x)

def rastrigin( x ):
  f = 0.0
  for c in x:
    f += rastrigin_1(c)
  return 10.0*len(x) + f

# Schwefel Function
def schwefel_1( x ):
	return -x * math.sin(math.sqrt(math.abs(x)))

def schwefel( x ):
  f = 0.0
  for c in x:
    f += schwefel(c)
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
  return -a*math.exp(-b*(p/D)**0.5) - math.exp(q/D)

##################### TEST PROBLEMS ####################
def RealTestProblem(f, D, EVALS, TRACE=False):
  if(f=='Rastrigin'): problem = PROBLEM('min', rastrigin, RealSpace(-5.12, 5.12, D), EVALS, TRACE)
  elif(f=='Schwefel'): problem = PROBLEM('min', schwefel, RealSpace(-500.0, 500.0, D), EVALS, TRACE)
  elif(f=='Griewank'): problem = PROBLEM('min', griewank, RealSpace(-600.0, 600.0, D), EVALS, TRACE)
  elif(f=='Rosenbrock'): problem = PROBLEM('min', rosenbrock_saddle, RealSpace(-2.048, 2.048, D), EVALS, TRACE)
  elif(f=='Ackley'): problem = PROBLEM('min', ackley, RealSpace(-32.768, 32.768, D), EVALS, TRACE)
  elif(f=='Sphere'): problem = PROBLEM('min', sphere, RealSpace(-5.12, 5.12, D), EVALS, TRACE)
  else: problem = PROBLEM('min', sphere, RealSpace(-5.12, 5.12, D), EVALS, TRACE)
  problem['optimum'] = 0.0
  return problem

##################### TEST PROBLEMS AS BINARY ####################
def Real2BinaryTestProblem(f, D, EVALS, BITSIZE = 32, TRACE=False):
  if(f=='Rastrigin'): problem = Real2BinaryPROBLEM('min', rastrigin, Real2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  elif(f=='Schwefel'): problem = Real2BinaryPROBLEM('min', schwefel, Real2Binary(-500.0, 500.0, D, BITSIZE), EVALS, TRACE)
  elif(f=='Griewank'): problem = Real2BinaryPROBLEM('min', griewank, Real2Binary(-600.0, 600.0, D, BITSIZE), EVALS, TRACE)
  elif(f=='Rosenbrock'): problem = Real2BinaryPROBLEM('min', rosenbrock_saddle, Real2Binary(-2.048, 2.048, D, BITSIZE), EVALS, TRACE)
  elif(f=='Ackley'): problem = Real2BinaryPROBLEM('min', ackley, Real2Binary(-32.768, 32.768, D, BITSIZE), EVALS, TRACE)
  elif(f=='Sphere'): problem = Real2BinaryPROBLEM('min', sphere, Real2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  else: problem = Real2BinaryPROBLEM('min', sphere, Real2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  problem['optimum'] = 0.0
  return problem