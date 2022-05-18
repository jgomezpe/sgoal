# The "GABO: Gene Analysis Base Optimization" algorithm 
# proposed by Professors Jonatan Gomez and Elizabeth Leon from 
# Universidad Nacional de Colombia 
# published in proceedings of the IEEE World Congress on Computational
# Intelligence - WCCI 2022
import random as rand
from sgoal import permutation
from sgoal import randbool
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

# Checks all the information about the gene's contribution to get the higher one,
# determines the best allele (0, 1, or None) for the gene  and the first time (in checking trials) it was reached
def best_allele(k):
  global contribution
  ba, bc = -1,-1
  for c in contribution[k]:
    if(c>0): oa,oc = 1,c
    elif(c<0): oa,oc =  0,-c
    else: oa,oc = 0 if randbool() else 1, 0
    if(oc > bc): ba, bc = oa, oc
  return ba

# Best by gene contribution
def genome_best_gene_contribution(f, evals, x, fx):
  if(evals>0):
    y = [best_allele(k) for k in range(len(x))] 
    fy = evaluate(f,y)
    x, y, fx, fy = pick(x, y, fx, fy)
    evals -= 1
  return x, fx, evals
  
# initializes global variables
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

  xc = complement(x)
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
    yc = complement(y)
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

def init(D):
  global intron, coding, separable, contribution
  separable = [True for k in range(D)]
  contribution = [[] for k in range(D)]
  intron = [k for k in range(D)]
  coding = []

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