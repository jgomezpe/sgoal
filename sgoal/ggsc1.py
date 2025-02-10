# 1. A Generalization of The Global Search Algorithm for GA-easy functions 
# propossed by Das and Whitley 1991, which tries only order 1-schemas, see
# R. Das and L. D. Whitley, “The only challenging problems are deceptive: 
# Global search by solving order-1 hyperplanes,
# in Proceedings of the 4th International Conference on Genetic Algorithms, San Diego, CA, USA,
# July 1991 (R. K. Belew and L. B. Booker, eds.), pp. 166– 173, Morgan Kaufmann, 1991
# 2. A Generalization of The Global Search Algorithm with complement propossed by 
# G. Venturini 1995  in “Ga consistently deceptive functions are not challenging problems” in 
# First International Conference on Genetic Algorithms in Engineering Systems: 
# Innovations and Applications, pp. 357–364, 1995.
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
from sgoal.binary import complement


# A Generalization of The Global Search Algorithm for GA-easy functions 
# propossed by Das and Whitley 1991, which tries only order 1-schemas, see
# R. Das and L. D. Whitley, “The only challenging problems are deceptive: 
# Global search by solving order-1 hyperplanes,
# in Proceedings of the 4th International Conference on Genetic Algorithms, San Diego, CA, USA,
# July 1991 (R. K. Belew and L. B. Booker, eds.), pp. 166– 173, Morgan Kaufmann, 1991.
# Our generalization allows to check a good 1-schema after some CHECK evaluations
# instead of checing it at the end of the allowed number of fitness evaluations
# f: Function to be optimized
def nextPopGGS1(P, fP, sgoal):
  if(sgoal.check == -1): 
    sgoal.check = sgoal.evals
  if((sgoal.count+sgoal.delta)%sgoal.check != 0):
    y = sgoal.space.get(1)
    return y, sgoal.evalone(y)

# Computes schemata information
  M = len(sgoal.S1)
  D = sgoal.space.D
  C = [[0 for k in range(D)], [0 for k in range(D)]]
  fH = [[0 for k in range(D)], [0 for k in range(D)]]
  for k in range(D):
    for i in range(M):
      fH[sgoal.S1[i][k]][k] += sgoal.fS1[i]
      C[sgoal.S1[i][k]][k] += 1
  # Generates a candidate solution with the best genes
  y = []
  for k in range(D):
    if( sgoal.minimize ):
      y.append( 1 if(fH[1][k]/C[1][k] < fH[0][k]/C[0][k]) else 0 )
    else:
      y.append( 1 if(fH[1][k]/C[1][k] > fH[0][k]/C[0][k]) else 0 )
  fy = SGoal.evalone(sgoal, y)
  return y, fy

# Chunks Global Search Algorithm. Computes GGS1 every CHECK function evaluations
# problem: Problem to solve
# CHECK: Function evaluations for computing 'the best so far' solution
class GGS1(SGoal):
  def __init__(self, problem, CHECK):
    SGoal.__init__(self, problem, nextPopGGS1)
    self.S1 = []
    self.fS1 = []
    self.N = 1
    self.check = CHECK
    self.delta = 1

  # Evals one candidate solution
  def evalone(self, x):
    fx = SGoal.evalone(self, x)
    self.S1.append(x)
    self.fS1.append(fx)
    return fx


# Global Search Algorithm 
# problem: Problem to solve
def GS1(problem): return GGS1(problem, -1)

# Chunks Global Search Algorithm. Computes GS1 every D function evaluations
# problem: Problem to solve
def DGS1(problem): return GGS1(problem, problem['space'].D)


# A Generalization of The Global Search Algorithm with complement propossed by 
# G. Venturini 1995
# Applies the GS1 and compares the obtained solution with its complement and the
# best candidate solution of the S1 set (same as the best solution found), see
# G. Venturini, “Ga consistently deceptive functions are not challenging problems”
# in First International Conference on Genetic Algorithms in Engineering Systems: 
# Innovations and Applications, pp. 357–364, 1995.
# Our generalization checks the complement after some fitness evaluations
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Initial point, a random point. It is required for defining the dimension of the space
# fx: f value at point x (if provided)
def nextPopGGSC1(P, fP, sgoal):
    y, fy = nextPopGGS1(P, fP, sgoal)
    if((sgoal.count+1)%sgoal.check == 0):
      y = complement(sgoal.result['x'])
      fy = SGoal.evalone(sgoal, y)
    return y, fy

# Chunks Global Search Algorithm with complement. Computes GGSC1 every CHECK function evaluations
# problem: Problem to solve
# CHECK: Function evaluations for computing 'the best so far' solution
def GGSC1(problem, CHECK):
    sgoal = GGS1(problem, CHECK)
    sgoal.delta = 2
    return sgoal

# Global Search Algorithm with complement. Computes GGSC1 every D function evaluations
# problem: Problem to solve
def GSC1(problem): return GGSC1(problem, -1)

# Chunks Global Search Algorithm with COmplement. Computes GS1 every D function evaluations
# problem: Problem to solve
def DGSC1(problem): return GGSC1(problem, problem['space'].D)
