# The "GABO: Gene Analysis Base Optimization" algorithm as proposed by
# J. Gomez and E. Leon, "Gabo: Gene Analysis Bitstring Optimization," 
# 2022 IEEE Congress on Evolutionary Computation (CEC), Padua, Italy, 2022, 
# pp. 1-8, doi: 10.1109/CEC55065.2022.9870237.
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
import random as rand
from sgoal.core import caneval
from sgoal.core import SPSGoal
from sgoal.util import permutation
from sgoal.binary import multiflip
from sgoal.binary import complement
from sgoal.binary import flip
from sgoal.core import initPop
from sgoal.core import init


# Computes contribution information (relative to a value 1), i.e., some change 
# in the f value. Minimization version
# x: A candidate solution
# y: The candidate solution with the k-th bit flipped 
# fx: f value of x
# fy: f value of y
# k: Gene's locus
def min_C(x, fx, fy, k, sgoal):
  c = fy-fx
  if(x[k]==0): c = -c 
  sgoal['C'][k].append(c)
  return c

# Computes contribution information (relative to a value 1), i.e., some change 
# in the f value. Maximization version
# x: A candidate solution
# y: The candidate solution with the k-th bit flipped 
# fx: f value of x
# fy: f value of y
# k: Gene's locus
def max_C(x, fx, fy, k, sgoal):
  c = fx-fy
  if(x[k]==0): c = -c 
  sgoal['C'][k].append(c)
  return c

# Evals contribution for each gene and gets the best according to improvements
def allelesCheck(x, fx, sgoal):
  pick, flip = sgoal['pick'], sgoal['flip']
  if(sgoal['minimize']): C = min_C
  else: C = max_C

  P = permutation(len(x))
  for k in P:
    if(not caneval(sgoal)): return x, fx
    y, fy = flip(x,fx,k)
    x, fx, y, fy = pick(x, fx, y, fy)
    C(x, fx, fy, k, sgoal)
  return x, fx

# Checks if a gene looks like intron
def intronLike(C):
  for k in C:
    if(k!=0): return False
  return True

# Splits genes into intron like and coding like genes according to computed contributions
def split(C):
  intron = []
  coding = []
  for k in range(len(C)):
    if(intronLike(C[k])): intron.append(k)
    else: coding.append(k)
  return intron, coding

# Initial Intron/Coding genes splitting
def ICSplit(x, fx, sgoal):
  x, fx = allelesCheck(x, fx, sgoal)
  sgoal['intron'], sgoal['coding'] = split(sgoal['C'])
  return x, fx

# Intron Only Search Algorithm
def IOSA(x, fx, sgoal):
  pick, flip, intron, nonintron = sgoal['pick'], sgoal['flip'], sgoal['intron'], sgoal['nonintron']
  if(sgoal['minimize']): C = min_C
  else: C = max_C
  N = len(intron)
  rand.shuffle(intron)
  i=0
  while(i<N):
    if(not caneval(sgoal)): return x, fx
    k = intron[i] # Picks and analyzes one intron-like locus
    y, fy = flip(x, fx, k)
    x, fx, y, fy = pick(x, fx, y, fy)
    if( C(x, fx, fy, k, sgoal) != 0 ):
      nonintron.append(k)
      intron.pop(i)
      N-=1
    else:
      i += 1
  return x, fx

# Checks all the information about the gene's contribution to get the higher one,
# determines the best allele (0, 1, or None) for the gene
# k: Locus (gene's position)
# b: Current value of the gene (allele)
def bestAllele(k, b, C):
  a, c = b, 0
  for ci in C[k]:
    if(ci>c): a, c = 1, ci
    elif(-ci>c): a, c =  0, -ci
  return a

# Generates a candidate solution with the alleles in the value having the highest contribution
def bestAlleles(x, sgoal):
  C = sgoal['C']
  return [bestAllele(k, x[k], C) for k in range(len(x))]

# Creates a genome with its best genes according to gene contribution and return the best between it and the original genome
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: current geneme, used for completing intron-like genes
# fx: Fitness value of the current genome
def bestByContribution(x, fx, sgoal):
  #if(not self.stop()):
  f, pick = sgoal['f'], sgoal['pick']
  y = bestAlleles(x, sgoal)
  fy = f(y)
  x, fx, y, fy = pick(x, fx, y, fy)
  return x, fx
  
# Coding Only Search Algorithm
def COSA( x, fx, sgoal ):
  pick, flip, multiflip, coding = sgoal['pick'], sgoal['flip'], sgoal['multiflip'], sgoal['coding']
  if(sgoal['minimize']): C = min_C
  else: C = max_C

  if(len(coding)==0 or not caneval(sgoal)): return x, fx

  N = len(coding)
  xc, fxc = multiflip(x, fx, coding)
  x, fx, xc, fxc = pick(x, fx, xc, fxc)

  # Considers locus by locus in a random fashion
  perm = permutation(N)
  for i in perm:
    k = coding[i]
    if(not caneval(sgoal)): return x, fx
    y, fy = flip(x, fx, k)
    yc, fyc = multiflip(y, fy, coding)
    C(x, fx, fy, k, sgoal)
    C(xc, fxc, fyc, k, sgoal)
    w = x
    y, fy, yc, fyc = pick( y, fy, yc, fyc )
    x, fx, y, fy = pick( x, fx, y, fy )
    if(w!=x): xc, fxc = yc, fyc
  x, fx = bestByContribution(x, fx, sgoal)
  return x, fx

  
############ GABO: Gene Analysis Bitstring Optimization #############
# Next Population method of the GABO Algorithm
def next( x, fx, sgoal ):
  x, fx = IOSA(x, fx, sgoal)
  x, fx = COSA(x, fx, sgoal)
  return x, fx

# Init individual method of the GABO Algorithm
def initGABO(sgoal):
  x, fx = init(sgoal)
  x, fx = ICSplit(x, fx, sgoal)
  x, fx = bestByContribution(x, fx, sgoal)
  x, fx = COSA(x, fx, sgoal)
  return x, fx

# Flip a bit -> variation form
def sflip(x, fx, k, sgoal): 
  y = flip(x, k)
  fy = sgoal['f'](y)
  return y, fy

# Multi Flip bits -> variation form
def mflip(x, fx, k, sgoal): 
  y = multiflip(x, k)
  fy = sgoal['f'](y)
  return y, fy

# GABO Algorithm Configuration. Extends a Binary Problem with the follwoing keys
#   'C': Contribution information (A vector with all alleles computed contributions) 
#   'intron': Array with intron like allele indices
#   'nonintron': Array with non-intron like allele indices according to IOSA
#   'coding': Array with coding allele indices
#   'flip': Bit flip method --> variation form (by default sets sflip)
#   'flip': MultiBit flip method --> variation form (by default sets mflip)
def GABOConfig(problem):
  D = problem['D']
  problem['C'] = [[] for k in range(D)]
  problem['intron'] = [k for k in range(D)]
  problem['nonintron'] = []
  problem['coding'] = []
  if('flip' not in problem): problem['flip'] = lambda x, fx, k: sflip(x, fx, k, problem)
  if('multiflip' not in problem): problem['multiflip'] = lambda x, fx, k: mflip(x, fx, k, problem)
  return SPSGoal(problem)

# GABO Algorithm.
def GABO(problem):
  sgoal = GABOConfig(problem)
  sgoal['init'] = lambda : initGABO(sgoal)
  sgoal['next'] = lambda x, fx : next(x, fx, sgoal)
  return sgoal

################### ACIA: Alleles Contribution Initialization Algorithm (for any Binary SGOAL) using GABO ideas ######################
# Init individual method of the GABO Algorithm
def nextACIA( x, fx, sgoal ):
  f, pick = sgoal['f'], sgoal['pick']
  x, fx = allelesCheck(x, fx, sgoal)
  y, fy = bestByContribution(x, fx, sgoal)
  yc = complement(y)
  fyc = f(yc)
  yc, fyc = allelesCheck(yc, fyc, sgoal)
  y, fy, yc, fyc = pick(y, fy, yc, fyc)
  x, fx, y, fy = pick(x, fx, y, fy)
  y, fy = bestByContribution(x, fx, sgoal)
  x, fx, y, fy = pick(x, fx, y, fy)
  return x, fx

# ACIA
def ACIA(problem):
  D = problem['D']
  problem['EVALS'] = 2*D + 4
  problem['C'] = [[] for k in range(D)]
  if('flip' not in problem): problem['flip'] = lambda x, fx, k: sflip(x, fx, k, problem)
  if('multiflip' not in problem): problem['multiflip'] = lambda x, fx, k: mflip(x, fx, k, problem)
  problem['next'] = lambda x, fx: nextACIA(x, fx, problem)
  return SPSGoal(problem)
