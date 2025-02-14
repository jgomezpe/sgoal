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
from sgoal.core import SGoal
from sgoal.core import permutation
from sgoal.binary import flip
from sgoal.binary import multiflip
from sgoal.binary import complement
from sgoal.core import basicInitPop
from sgoal.core import basicStop

# Init GABO attributes
def initGABOInfo(sgoal):
  sgoal.contribution = []
  D = sgoal.space.D
  sgoal.contribution = [[] for k in range(D)]
  sgoal.candidate = [[] for k in range(D)]
  sgoal.intron = [k for k in range(D)]
  sgoal.nonintron = []
  sgoal.coding = []

# Computes contribution information (relative to a value 1), i.e., some change 
# in the f value
# x: A candidate solution
# y: The candidate solution with the k-th bit flipped 
# fx: f value of x
# fy: f value of y
# k: Gene's locus
def C(x, fx, fy, k, sgoal):
  if( sgoal.minimize ): c = fy-fx
  else: c = fx-fy
  if(x[k]==0): c = -c 
  sgoal.contribution[k].append(c)
  sgoal.candidate[k].append(x)
  return c
  
# Intron/Coding genes splitting
def ICSplit(x, fx, sgoal):
  P = permutation(len(x))
  for k in P:
    if(not sgoal.caneval()): return x, fx
    y = flip(x,k)
    fy = sgoal.evalone(y)
    x, fx, y, fy = sgoal.pick(x, fx, y, fy)
    if( C(x, fx, fy, k, sgoal) != 0 ): 
      sgoal.intron.remove(k)
      sgoal.coding.append(k)  
  return x, fx

# Intron Only Search Algorithm
def IOSA(x, fx, sgoal):
  N = len(sgoal.intron)
  rand.shuffle(sgoal.intron)
  i=0
  while(i<N):
    if(not sgoal.caneval()): return x, fx
    k = sgoal.intron[i] # Picks and analyzes one intron-like locus
    y = flip(x,k)
    fy = sgoal.evalone(y)
    x, fx, y, fy = sgoal.pick(x, fx, y, fy)
    if( C(x, fx, fy, k, sgoal) != 0 ):
      sgoal.nonintron.append(k)
      sgoal.intron.pop(i)
      N-=1
    else:
      i += 1
  return x, fx

# Checks if an allele looks like separable by considering its computed contributions
def checkseparable(C):
  c = C[0]
  for i in range(1,len(C)):
    if(c!=C[i]):
      return False
  return True

# Splits coding alleles in separable/nonseparable according to their contributions
def analize(sgoal):
  separable = []
  nonseparable = []
  for k in sgoal.coding:
    if(checkseparable(sgoal.contribution[k])):
      separable.append(k)
    else:
      nonseparable.append(k)
  for k in sgoal.coding:
    b = []
    for l in sgoal.coding:
      if(l!=k and sgoal.candidate[k][0][l]==sgoal.candidate[k][1][l]):
        b.append(l)
    x = sgoal.candidate[k][0]
    y = flip(x, k)
    xc = multiflip(x, b)
    fxc = sgoal.evalone(xc)
    yc = multiflip(y, b)
    fyc = sgoal.evalone(yc)
    C(xc, fxc, fyc, k, sgoal)
  separable = []
  nonseparable = []
  for k in sgoal.coding:
    if(checkseparable(sgoal.contribution[k])):
      separable.append(k)
    else:
      nonseparable.append(k)
  return separable, nonseparable

# Checks all the information about the gene's contribution to get the higher one,
# determines the best allele (0, 1, or None) for the gene
# k: Locus (gene's position)
# b: Current value of the gene (allele)
def best_allele(k, b, sgoal):
  a, c = b, 0
  for ci in sgoal.contribution[k]:
    if(ci>c): a, c = 1, ci
    elif(-ci>c): a, c =  0, -ci
  return a

# Creates a genome with its best genes according to gene contribution and return the best between it and the original genome
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: current geneme, used for completing intron-like genes
# fx: Fitness value of the current genome
def genome_best_gene_contribution(x, fx, sgoal):
  #if(not self.stop()):
  y = [best_allele(k, x[k], sgoal) for k in range(len(x))] 
  fy = sgoal.evalone(y)
  x, fx, y, fy = sgoal.pick(x, fx, y, fy)
  return x, fx

# Computes the coding complement of a genome. 
# The coding complement flips each gene that is considered coding and maintains other genes the same
# x: Genome to be complemented in the coding-like genes
def coding_complement(x, sgoal):
  N = len(sgoal.coding)
  if(N>0):
    y = x.copy()
    for i in sgoal.coding:
      y[i] = 1 if y[i]==0 else 0
  else:
    y = complement(x)
  return y
  
# Coding Only Search Algorithm
def COSA( x, fx, sgoal ):
  if(len(sgoal.coding)==0 or not sgoal.caneval()): 
    return x, fx
  N = len(sgoal.coding)
  xc = coding_complement(x, sgoal)
  fxc = sgoal.evalone(xc)
  x, fx, xc, fxc = sgoal.pick(x, fx, xc, fxc)

  # Considers locus by locus in a random fashion
  perm = permutation(N)
  for i in perm:
    k = sgoal.coding[i]
    if(not sgoal.caneval(2)): 
      return x, fx
    y = flip(x,k)
    fy = sgoal.evalone(y)
    yc = coding_complement(y, sgoal)
    fyc = sgoal.evalone(yc)    
    cx = C(x, fx, fy, k, sgoal)
    cxc = C(xc, fxc, fyc, k, sgoal)
    #if(cx!=cxc):
    #  self.nonseparable[k] = 1
    w = x
    y, fy, yc, fyc = sgoal.pick( y, fy, yc, fyc )
    x, fx, y, fy = sgoal.pick( x, fx, y, fy )
    if(w!=x): xc, fxc = yc, fyc
  x, fx = genome_best_gene_contribution(x, fx, sgoal)
  return x, fx

  
############ GABO: Gene Analysis Bitstring Optimization #############
# Next Population method of the GABO Algorithm
def nextPopGABO( P, fP, sgoal ):
  x, fx = IOSA(P, fP, sgoal)
  x, fx = COSA(x, fx, sgoal)
  return x, fx

# Init Population method of the GABO Algorithm
def initPopGABO(sgoal):
  initGABOInfo(sgoal)
  x, fx = basicInitPop(sgoal)
  x, fx = ICSplit(x, fx, sgoal)
  x, fx = genome_best_gene_contribution(x, fx, sgoal)
  x, fx = COSA(x, fx, sgoal)
  return x, fx

# GABO Algorithm
# problem: BitArray problem to solve
def GABO(problem):
  return SGoal(problem, nextPopGABO, initPopGABO, basicStop)

################### GABO2 ######################
# Analyzes coding alleles
def codingcheck(x, fx, sgoal):
  P = permutation(len(sgoal.coding))
  for i in P:
    if(not sgoal.caneval()): return x, fx
    k = sgoal.coding[i]
    y = flip(x,k)
    fy = sgoal.evalone(y)
    x, fx, y, fy = sgoal.pick(x, fx, y, fy)
    C(x, fx, fy, k, sgoal)
  return x, fx

# Init population methof of the GABO2 Algorithm
def initPopGABO2(sgoal):
  initGABOInfo(sgoal)
  x, fx = basicInitPop(sgoal)
  x, fx = ICSplit(x, fx, sgoal)
  y, fy = genome_best_gene_contribution(x, fx, sgoal)
  yc = coding_complement(y, sgoal)
  fyc = sgoal.evalone(yc)
  yc, fyc = codingcheck(yc, fyc, sgoal)
  y, fy, yc, fyc = sgoal.pick(y, fy, yc, fyc)
  x, fx, y, fy = sgoal.pick(x, fx, y, fy)
  y, fy = genome_best_gene_contribution(x, fx, sgoal)
  x, fx, y, fy = sgoal.pick(x, fx, y, fy)
  sgoal.coding.sort()
  sgoal.intron.sort()
  sep, nonsep = analize(sgoal)
  print('sep:', sep)
  print('nonsep:', nonsep)
  for i in range(len(sgoal.candidate[7])):
    for k in sgoal.candidate[7][i]:
      print(k, sep='', end='')
    print(' ', sgoal.contribution[7][i])
  return x, fx

def GABO2Stop(sgoal):
  return basicStop(sgoal) or len(sgoal.intron)==0

# GABO2 Algorithm
# problem: Problem to solve
def GABO2(problem):
  return SGoal(problem, IOSA, initPopGABO2, GABO2Stop)
