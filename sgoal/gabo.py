# Copyright (c)
# Authors: Jonatan Gomez and Elizabeth León  
# E-mails: jgomezpe@unal.edu.co and eleonguz@unal.edu.co
# All rights reserved.
# The "GABO: Gene Analysis Base Optimization" algorithm 
# proposed by Professors Jonatan Gomez and Elizabeth Leon from 
# Universidad Nacional de Colombia 
# published in proceedings of the IEEE World Congress on Computational
# Intelligence - WCCI 2022
import random as rand
from sgoal import SGoal
from sgoal import permutation
from binary import flip
from binary import complement

  
############ GABO: Gene Analysis Bitstring Optimization #############
# Gene information
class GABO(SGoal):
  def __init__(self, problem):
    SGoal.__init__(self, problem)
    self.contribution = []
    self.intron = []
    self.nonintron = []
    self.separable = []
    self.coding = []

  # Initializes the gene information 
  # D: Bitstring's length (number of genes)
  def init(self, MAXEVALS, TRACE):
    D = self.space.D
    x, fx = SGoal.init(self, MAXEVALS, TRACE)
    self.nonseparable = {}
    self.contribution = [[] for k in range(D)]
    self.intron = [k for k in range(D)]
    return x, fx 

  # Computes contribution information (relative to a value 1), i.e., some change 
  # in the f value
  # x: A candidate solution
  # y: The candidate solution with the k-th bit flipped 
  # fx: f value of x
  # fy: f value of y
  # k: Gene's locus
  def C(self, x, fx, fy, k):
    if( self.minimize ): c = fy-fx
    else: c = fx-fy
    if(x[k]==0): c = -c 
    self.contribution[k].append(c)
    return c

  # Checks all the information about the gene's contribution to get the higher one,
  # determines the best allele (0, 1, or None) for the gene
  # k: Locus (gene's position)
  # b: Current value of the gene (allele)
  def best_allele(self, k, b):
    a, c = b, 0
    for ci in self.contribution[k]:
      if(ci>c): a, c = 1, ci
      elif(-ci>c): a, c =  0, -ci
    return a

  # Creates a genome with its best genes according to gene contribution and return the best between it and the original genome
  # f: Function to be optimized
  # evals: Maximum number of fitness evaluations
  # x: current geneme, used for completing intron-like genes
  # fx: Fitness value of the current genome
  def genome_best_gene_contribution(self, x, fx):
    if(not self.stop()):
      y = [self.best_allele(k, x[k]) for k in range(len(x))] 
      fy = self.evalone(y)
      x, fx, y, fy = self.pick(x, fx, y, fy)
    return x, fx


  # Initializes global variables and does an initial gene analysis
  # f: Function to be optimized
  # evals: Maximum number of fitness evaluations
  # x: Current geneme
  # fx: Fitness value of the current genome
  def ICSplit(self, x, fx):
    P = permutation(len(x))
    for k in P:
      if(not self.caneval()): return x, fx
      y = flip(x,k)
      fy = self.evalone(y)
      x, fx, y, fy = self.pick(x, fx, y, fy)
      if( self.C(x, fx, fy, k) != 0 ): 
        self.intron.remove(k)
        self.coding.append(k)  
    return self.genome_best_gene_contribution(x, fx)
  
  def codingseparable(self):
    for k in self.coding:
          if(k in self.nonseparable): return False
    return True
  
  # Gabo's stop condition
  def stop(self):
    return SGoal.stop(self) or (self.codingseparable() and len(self.intron)==0)

  def IOSA( self, x, fx):
    N = len(self.intron)
    rand.shuffle(self.intron)
    i=0
    while(i<N):
      if(not self.caneval()): return x, fx
      k = self.intron[i] # Picks and analyzes one intron-like locus
      y = flip(x,k)
      fy = self.evalone(y)
      x, fx, y, fy = self.pick(x, fx, y, fy)
      if( self.C(x, fx, fy, k) != 0 ):
        self.nonintron.append(k)
        self.intron.pop(i)
        N-=1
      else:
        i += 1
    return x, fx


  # Computes the coding complement of a genome. 
  # The coding complement flips each gene that is considered coding and maintains other genes the same
  # x: Genome to be complemented in the coding-like genes
  def coding_complement(self, x):
    N = len(self.coding)
    if(N>0):
      y = x.copy()
      for i in self.coding:
        y[i] = 1 if y[i]==0 else 0
    else:
      y = complement(x)
    return y
    
  # Gene characterization algorithm
  # genome: An array with each gene information, see Gene class
  # f: Function to be optimized
  # x: initial point
  # fx: f value at point x (if provided)
  # evals: Maximum number of fitness evaluations
  def COSA( self, x, fx ):
    if(len(self.coding)==0 or not self.caneval()): 
      return x, fx
    N = len(self.coding)
    xc = self.coding_complement(x)
    fxc = self.evalone(xc)
    x, fx, xc, fxc = self.pick(x, fx, xc, fxc)

    # Considers locus by locus in a random fashion
    perm = permutation(N)
    for i in perm:
      k = self.coding[i]
      if(not self.caneval(2)): 
        return x, fx
      y = flip(x,k)
      fy = self.evalone(y)
      yc = self.coding_complement(y)
      fyc = self.evalone(yc)    
      cx = self.C(x, fx, fy, k)
      cxc = self.C(xc, fxc, fyc, k)
      if(cx!=cxc):
        self.nonseparable[k] = 1
      w = x
      y, fy, yc, fyc = self.pick( y, fy, yc, fyc )
      x, fx, y, fy = self.pick( x, fx, y, fy )
      if(w!=x): xc, fxc = yc, fyc
    x, fx = self.genome_best_gene_contribution(x, fx)
    return x, fx


  # f: Function to be optimized
  # x: initial point
  # evals: Maximum number of fitness evaluations
  # fx: f value at point x (if provided)
  def run( self, MAXEVALS, TRACE=False):
    x, fx = self.init(MAXEVALS, TRACE)
    x, fx = self.ICSplit(x, fx)
    x, fx = self.COSA(x, fx)
    while(not self.stop()):
      x, fx = self.IOSA(x, fx)
      x, fx = self.COSA(x, fx)
    self.tracing(x,fx)
    if('evals' not in self.result):
      self.result['evals'] = self.count
    return self.result
