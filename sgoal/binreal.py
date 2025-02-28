# Real problem to Binary (BitArray) problem definitions
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

from sgoal.core import PROBLEM
from sgoal.binary import Binary
from sgoal.real import rastrigin
from sgoal.real import sphere
from sgoal.real import rosenbrock_saddle
from sgoal.real import ackley
from sgoal.real import griewank
from sgoal.real import schwefel

##################### INTERVAL TO BINARY #####################
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

# Interval to Binary Space
def Interval2Binary(min, max, BITSIZE=32):
  D = BITSIZE
  space = Binary(D)
  length = max - min
  space['grow'] = lambda x: grow1(x, min, length, BITSIZE)
  return space

##################### HYPERECTANGLE/HYPERCUBE TO BINARY #####################
# Grows a binary representation to real vector representation
def grow(x, min, length, SIZE=16):
  return [grow1(x, min[i], length[i], SIZE, i) for i in range(len(x)//SIZE)]

# HyperRectangle to Binary Space
def HyperRectangle2Binary(min, max, BITSIZE=32):
  D = len(min)*BITSIZE
  space = Binary(D)
  length = [max[i]-min[i] for i in range(len(min))]
  space['grow'] = lambda x: grow(x, min, length, BITSIZE)
  return space

# HyperCube to Binary Space
def HyperCube2Binary(min, max, D=2, BITSIZE=32):
  if(D<2): D=2
  return HyperRectangle2Binary([min for i in range(D)], [max for i in range(D)], BITSIZE)

# Realvalued Problem to Binary problem
def Real2BinaryPROBLEM(type, f, space, EVALS, TRACE=False):
  return PROBLEM(type, lambda x: f(space['grow'](x)), space, EVALS, TRACE)

##################### Real valued test problems as Binary problems ####################
def Real2BinaryTestProblem(f, D, EVALS, BITSIZE = 32, TRACE=False):
  if(D<2): D = 2
  if(f=='Rastrigin'): problem = Real2BinaryPROBLEM('min', rastrigin, HyperCube2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  elif(f=='Schwefel'): problem = Real2BinaryPROBLEM('min', schwefel, HyperCube2Binary(-500.0, 500.0, D, BITSIZE), EVALS, TRACE)
  elif(f=='Griewank'): problem = Real2BinaryPROBLEM('min', griewank, HyperCube2Binary(-600.0, 600.0, D, BITSIZE), EVALS, TRACE)
  elif(f=='Rosenbrock'): problem = Real2BinaryPROBLEM('min', rosenbrock_saddle, HyperCube2Binary(-2.048, 2.048, D, BITSIZE), EVALS, TRACE)
  elif(f=='Ackley'): problem = Real2BinaryPROBLEM('min', ackley, HyperCube2Binary(-32.768, 32.768, D, BITSIZE), EVALS, TRACE)
  elif(f=='Sphere'): problem = Real2BinaryPROBLEM('min', sphere, HyperCube2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  else: problem = Real2BinaryPROBLEM('min', sphere, HyperCube2Binary(-5.12, 5.12, D, BITSIZE), EVALS, TRACE)
  problem['optimum'] = 0.0
  return problem