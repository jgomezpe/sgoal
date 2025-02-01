# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co and eleonguz@unal.edu.co
# All rights reserved.
# The "GABO: Gene Analysis Base Optimization" algorithm 
# proposed by Professors Jonatan Gomez and Elizabeth Leon from sgoal.
# Universidad Nacional de Colombia 
# published in proceedings of the IEEE World Congress on Computational
# Intelligence - WCCI 2022
from sgoal.gabo import GABO

import random as rand
from sgoal.core import SGoal
from sgoal.core import permutation
from sgoal.binary import flip
from sgoal.binary import multiflip

class GABO2(GABO):
  def __init__(self, problem):
    GABO.__init__(self, problem)
    self.relation = []


  # Initializes the gene information 
  # D: Bitstring's length (number of genes)
  def init(self, MAXEVALS, TRACE):
    D = self.space.D
    x, fx = GABO.init(self, MAXEVALS, TRACE)
    self.relation = [[] for i in range(D)]
    return x, fx
  
  def REL_BIS(self, x, fx, y, fy, bit1, groups, start, end):
    if(not self.caneval(2)): 
      return x, fx, -1
    bit2 = []
    for k in range(start,end):
      bit2 += groups[k]
    xc = multiflip(y,bit2)
    fxc = self.evalone(xc)
    yc = flip(xc,bit1)
    fyc = self.evalone(yc)
    if(abs(fx-fxc)!=abs(fx-fy)+abs(fx-fyc)):
      m = (start+end)//2
      if(m==start): 
        return x, fx, k
      x, fx, k = self.REL_BIS(x, fx, y, fy, bit1, groups, start, m)
      if(k<0):
        return self.REL_BIS(x, fx, y, fy, bit1, groups, m, end)
      else:
        return x, fx, k
    else:
      return x, fx, -1    
  
  def stop3(self):
    return self.evals<self.count+3
  
  def REL(self, x, fx, genes):
    M = len(genes)
    if(M==0):
      return x, fx
    groups = [[genes[M-1]]]
    genes.pop()
    M -= 1
    i=M-1
    while(i>=0 and not self.stop3()):
      end = len(groups)
      bit1 = genes[i]
      y = flip(x,bit1)
      fy = self.evalone(y)
      x, fx, k = self.REL_BIS(x, fx, y, fy, bit1, groups, 0, end)
      if(k>=0):
        groups[k].append(genes[i])
      else:
        groups.append([genes[i]])
      genes.pop()
      i -= 1
    for l in groups:
      for i in l:
        for j in l:
          if(i!=j):
            self.relation[i].append(j)
    return x, fx
  

  def separable_check(self, x, fx, xc, fxc, set, nonsep):
    N = len(set)
    if(N==0): return x, fx, xc, fxc, set, nonsep

    set2 = []

    # Considers locus by locus in a random fashion
    perm = permutation(N)
    for i in perm:
      k = set[i]
      if(not self.caneval(2)): return  x, fx, xc, fxc, set, nonsep
      y = flip(x,k)
      fy = self.evalone(y)
      yc = self.coding_complement(y)
      fyc = self.evalone(yc)    
      cx = self.C(x, fx, fy, k)
      cxc = self.C(xc, fxc, fyc, k)
      if(cx==cxc):
        set2.append(k)
      else:
        nonsep.append(k)
      w = x
      y, fy, yc, fyc = self.pick( y, fy, yc, fyc )
      x, fx, y, fy = self.pick( x, fx, y, fy )
      if(w!=x): xc, fxc = yc, fyc
    return x, fx, xc, fxc, set2, nonsep
  
  # Gene characterization algorithm
  # genome: An array with each gene information, see Gene class
  # f: Function to be optimized
  # x: initial point
  # fx: f value at point x (if provided)
  # evals: Maximum number of fitness evaluations
  def InitialSLSplit( self, x, fx ):
    if(not self.caneval()): return x, fx
    xc = self.coding_complement(x)
    fxc = self.evalone(xc)
    x, fx, xc, fxc = self.pick(x, fx, xc, fxc)
    x, fx, xc, fxc, self.separable, self.coding = self.separable_check(x, fx, xc, fxc, self.coding, [])
    return self.genome_best_gene_contribution(x, fx)

  
  # Gene characterization algorithm
  # genome: An array with each gene information, see Gene class
  # f: Function to be optimized
  # x: initial point
  # fx: f value at point x (if provided)
  # evals: Maximum number of fitness evaluations
  def SOSA( self, x, fx, set ):
    if(not self.caneval()): return x, fx, set

    xc = self.coding_complement(x)
    fxc = self.evalone(xc)
    x, fx, xc, fxc = self.pick(x, fx, xc, fxc)
    
    N = len(set)
    M=0
    while(M!=N and N>0):
      x, fx, xc, fxc, set2, ccod = self.separable_check(x, fx, xc, fxc, set, self.coding)    
      set = set2
      M = N
      N = len(set)
    x, fx = self.genome_best_gene_contribution(x, fx)
    return x, fx, set

  # Gabo's stop condition
  def stop(self):
    return SGoal.stop(self) or len(self.intron)==0

  # f: Function to be optimized
  # x: initial point
  # evals: Maximum number of fitness evaluations
  # fx: f value at point x (if provided)
  def run( self, MAXEVALS, TRACE=False):
    x, fx = self.init(MAXEVALS, TRACE)
    x, fx = self.ICSplit(x, fx)
    x, fx = self.InitialSLSplit(x, fx)
    #print("separable COSA:",separable)
    #print("intron COSA:",intron)
    x, fx, self.separable = self.SOSA(x, fx, self.separable)
    x, fx, self.intron = self.SOSA(x, fx, self.intron)
    x, fx = self.REL(x, fx, self.coding)
    while(not self.stop()):
      N = len(self.nonintron)
      x, fx = self.IOSA(x, fx)
      M = len(self.nonintron)
      if(M>0 and (N==M or len(self.intron)==0)):
        x, fx = self.REL(x, fx, self.nonintron)

    R = []
    g = [True for i in range(self.space.D)]
    for i in range(self.space.D):
      if(g[i]):
        self.relation[i].insert(0,i)
        R.append(self.relation[i])
        for k in self.relation[i]:
          g[k] = False
    self.result['groups'] = R
    self.tracing(x,fx)
    if('evals' not in self.result):
      self.result['evals'] = self.count
    return self.result
