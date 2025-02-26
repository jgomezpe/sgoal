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
from sgoal.core import VRSGoal
from sgoal.core import evalPop
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
def variationGS1(x, fx, sgoal):
  f, getN, minimize = sgoal['f'], sgoal['getN'], sgoal['minimize']
  
  N = sgoal['EVALS'] - 1
  S1 = getN(N)
  fS1 = evalPop(S1, f)
  S1.append(x)
  fS1.append(fx)

# Computes schemata information
  M = len(S1)
  D = sgoal['D']
  C = [[0 for k in range(D)], [0 for k in range(D)]]
  fH = [[0 for k in range(D)], [0 for k in range(D)]]
  for k in range(D):
    for i in range(M):
      fH[S1[i][k]][k] += fS1[i]
      C[S1[i][k]][k] += 1
  # Generates a candidate solution with the best genes
  y = []
  if( minimize ):
    for k in range(D):
      y.append( 1 if(fH[1][k]/C[1][k] < fH[0][k]/C[0][k]) else 0 )
  else:
    for k in range(D):
      y.append( 1 if(fH[1][k]/C[1][k] > fH[0][k]/C[0][k]) else 0 )
  fy = f(y)
  return y, fy

# Global Search Algorithm 
# problem: Problem to solve
def GS1(problem):
  problem['next'] = lambda x, fx: variationGS1(x, fx, problem)
  return VRSGoal(problem)


# The Global Search Algorithm with complement propossed by G. Venturini 1995
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
def variationGSC1(x, fx, sgoal):
  sgoal['EVALS']-=1
  y, fy = variationGS1(x, fx, sgoal)
  x, fx = sgoal['best']['x'], sgoal['best']['f']
  sgoal['EVALS']+=1
  y = complement(x)
  fy = sgoal['f'](y)
  return y, fy

# Global Search Algorithm with complement.
# problem: Problem to solve
def GSC1(problem):
  problem['next'] = lambda x, fx : variationGSC1(x, fx, problem)
  return GS1(problem)