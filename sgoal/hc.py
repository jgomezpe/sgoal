# Classical Hill Climbing Algorithm
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
from sgoal.binary import bitmutation
from sgoal.binary import singlebitmutation
from sgoal.real import gaussianMutation

############### BINARY HILL  CLIMBING ###############
# Classical Hill Climbing Algorithm with neutral mutation for BitArray problems. Uses bitmutation as variation operator
# problem: Problem to solve
def BinaryHC(problem): 
  if( 'variation' not in problem ): problem['variation'] = bitmutation 
  return VRSGoal(problem)

# The HC algorithm suggested by Richard Palmer, that Forrest and
# Mitchell named as "random mutation hill-climbing" (RMHC), see
# M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?”
# Santa Fe Institute, Working Papers, 01 1993.
# problem: Problem to solve
def RMHC(problem): 
  if('variation' not in problem): problem['variation'] = singlebitmutation 
  return VRSGoal(problem)

############### REAL HILL  CLIMBING ###############
# Classical Hill Climbing Algorithm for Real problems. Uses Gaussian mutation with sigma=0.2 as variation operator
# problem: Problem to solve
def RealHC(problem): 
  if( 'variation' not in problem ): problem['variation'] = gaussianMutation(problem)
  return VRSGoal(problem) 
