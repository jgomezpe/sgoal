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
from sgoal import permutation
from sgoal import pick
from sgoal import MAXIMIZE
from sgoal import evaluate
from binary import flip
from binary import complement

############ GABO: Gene Analysis Bitstring Optimization #############
# Gene information
contribution = []
intron = []
separable = []
coding = []

# Initializes the gene information 
# D: Bitstring's length (number of genes)
def init(D):
  global intron, coding, separable, contribution
  separable = [True for k in range(D)]
  contribution = [[] for k in range(D)]
  intron = [k for k in range(D)]
  coding = []

# Computes contribution information (relative to a value 1), i.e., some change 
# in the f value
# x: A candidate solution
# y: The candidate solution with the k-th bit flipped 
# fx: f value of x
# fy: f value of y
# k: Gene's locus
def C(x, y, fx, fy, k):
  global contribution
  if( MAXIMIZE ): c = fx-fy
  else: c = fy-fx
  if(x[k]==0): c = -c 
  contribution[k].append(c)
  return c

# Computes the coding complement of a genome. 
# The coding complement flips each gene that is considered coding and maintains other genes the same
# x: Genome to be complemented in the coding-like genes
def coding_complement(x):
  global coding
  N = len(coding)
  if(N>0):
    y = x.copy()
    for i in coding:
      y[i] = 1-y[i]
  else:
    y = complement(x)
  return y
      
# Checks all the information about the gene's contribution to get the higher one,
# determines the best allele (0, 1, or None) for the gene
# k: Locus (gene's position)
# b: Current value of the gene (allele)
def best_allele(k, b):
  global contribution
  a, c = b, 0
  for ci in contribution[k]:
    if(ci>c): a, c = 1, ci
    elif(-ci>c): a, c =  0, -ci
  return a

# Creates a genome with its best genes according to gene contribution and return the best between it and the original genome
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: current geneme, used for completing intron-like genes
# fx: Fitness value of the current genome
def genome_best_gene_contribution(f, evals, x, fx):
  if(evals>0):
    y = [best_allele(k, x[k]) for k in range(len(x))] 
    fy = evaluate(f,y)
    x, y, fx, fy = pick(x, y, fx, fy)
    evals -= 1
  return x, fx, evals
  
# Initializes global variables and does an initial gene analysis
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: Current geneme
# fx: Fitness value of the current genome
def ICSplit(f, evals, x, fx):
  global intron, coding, separable, contribution
  P = permutation(len(x))
  for k in P:
    if(evals==0): return x, fx, evals
    y = flip(x,k)
    fy = evaluate(f,y)
    evals -= 1    
    x, y, fx, fy = pick(x, y, fx, fy)
    if( C(x, y, fx, fy, k) != 0 ): 
      intron.remove(k)
      coding.append(k)  
  return genome_best_gene_contribution(f, evals, x, fx)

# Introns only search algorithm
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# x: initial point
# fx: f value in the x point
def IOSA( f, evals, x, fx):
  global intron
  N = len(intron)
  for i in range(N):
    if(evals==0): return x, fx, evals
    a = rand.randint(0,N-1)
    k = intron[a] # Picks and analyzes one intron-like locus
    y = flip(x,k)
    fy = evaluate(f,y)
    evals -= 1    
    x, y, fx, fy = pick(x, y, fx, fy)
    if( C(x, y, fx, fy, k) != 0 ):
      intron.pop(a)
      N-=1
  return x, fx, evals

# Gene characterization algorithm
# genome: An array with each gene information, see Gene class
# f: Function to be optimized
# x: initial point
# fx: f value at point x (if provided)
# evals: Maximum number of fitness evaluations
def COSA( f, evals, x, fx ):
  global coding, separable
  N = len(coding)
  if(N==0 or evals==0): return x, fx, evals

  xc = coding_complement(x)
  fxc = evaluate(f,xc)
  evals -= 1
  x, xc, fx, fxc = pick(x, xc, fx, fxc)

  # Considers locus by locus in a random fashion
  perm = permutation(N)
  for i in perm:
    k = coding[i]
    if(evals<2): return x, fx, evals
    y = flip(x,k)
    fy = evaluate(f,y)
    yc = coding_complement(y)
    fyc = evaluate(f,yc)    
    evals -= 2
    cx = C(x, y, fx, fy, k)
    cxc = C(xc, yc, fxc, fyc, k)
    separable[k] = separable[k] and (cx==cxc)
    w = x
    y, yc, fy, fyc = pick( y, yc, fy, fyc )
    x, y, fx, fy = pick( x, y, fx, fy )
    if(w!=x): xc, fxc = yc, fyc
  x, fx, evals = genome_best_gene_contribution(f, evals, x, fx)
  return x, fx, evals

# Gabo's stop condition
def stop():
  global coding, intron, separable
  for k in coding:
    if(not separable[k]): return False
  return len(intron)==0
  
# f: Function to be optimized
# x: initial point
# evals: Maximum number of fitness evaluations
# fx: f value at point x (if provided)
def GABO( f, evals, x, fx=None):
  init(len(x))
  if(evals>0 and not fx):
    fx = evaluate(f,x)
    evals -= 1
  x, fx, evals = ICSplit(f, evals, x, fx)
  flag = True
  while(evals>0 and flag):
    x, fx, evals = IOSA(f, evals, x, fx)
    x, fx, evals = COSA(f, evals, x, fx)
    flag = not stop()
  return x, fx 