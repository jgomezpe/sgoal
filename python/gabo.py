import random as rand
from sgoal import  PRINT
from sgoal import permutation
from sgoal import randbool
from sgoal import pick
from sgoal import MAXIMIZE
from sgoal import evaluate
from binary import bitstring
from binary import flip
from binary import complement
from binary import for_all

############ GABO: Gene Analysis Bitstring Optimization #############
# Gene information
contribution = []
intron = []
separable = []

# initializes global variables
def init_gene(D):
  global contribution, intron, separable
  contribution = [[] for i in range(D)]
  intron = [True for i in range(D)]
  separable = [True for i in range(D)]

# Computes contribution information (relative to a value 1), i.e., some change 
# in the f value
# x: A candidate solution
# y: The candidate solution with the k-th bit flipped 
# fx: f value of x
# fy: f value of y
# k: Gene's locus
def C(x, y, fx, fy, k):
  global intron, contribution
  if( MAXIMIZE ): c = fx-fy
  else: c = fy-fx
  intron[k] = intron[k] and c==0
  if(x[k]==0): c = -c 
  contribution[k].append(c)
  return c

# Gene characterization algorithm
# genome: An array with each gene information, see Gene class
# f: Function to be optimized
# x: initial point
# fx: f value at point x (if provided)
# evals: Maximum number of fitness evaluations
def GCA( f, evals, x, fx=None ):
  global separable
  D = len(x) # Space dimension
  if(evals>0 and not fx):
    fx = evaluate(f, x)
    evals -=1
  if( evals>0 ):
    # The complement candidate solution (used for determining locus separability)
    xc = complement(x)
    fxc = evaluate(f,xc)
    x, xc, fx, fxc = pick(x, xc, fx, fxc)
    evals -= 1

  # Considers locus by locus (shuffles loci for reducing order effect)
  perm = permutation(D)
  a=0
  while(evals>=2 and a<D):
    k = perm[a]
    y = flip(x,k)
    fy = evaluate(f,y)
    yc = complement(y)
    fyc = evaluate(f,yc)    
    cx = C(x, y, fx, fy, k)
    cxc = C(xc, yc, fxc, fyc, k)
    separable[k] = separable[k] and (cx==cxc)
    w = x
    y, yc, fy, fyc = pick( y, yc, fy, fyc )
    x, y, fx, fy = pick( x, y, fx, fy )
    if(w!=x): xc, fxc = yc, fyc
    a += 1
    evals -= 2
  if(PRINT): print('Best * GCA trial:', fx, x)
  print('Evals:',evals)
  return x, fx, evals

# Checks the gene's contribution information to determine if the best bit value (allele)
# for the gene is 0: c<0, or 1: c>0. If c=0 then a gene behaves like an intron (neutral) so 
# bit is set to a random value. Returns a non-negative contribution (c) accordingly
def allele(c):
  if(c>0): return 1,c
  if(c<0): return 0,-c
  bit = 0 if randbool() else 1
  return bit, c

# Checks all the information about the gene's contribution to get the higher one,
# determines the best allele (0, 1, or None) for the gene  and the first time (in checking trials) it was reached
def best_allele(k):
  global contribution
  alle, cont = allele(contribution[k][0])
  trial = 0
  for i in range(1,len(contribution[k])):
    a, c = allele(contribution[k][i])
    if(c > cont): alle, cont, trial = a, c, i
  return alle, cont, trial

# Creates a candidate solution using locus contributions.
# For each locust gets the value of the bit (0, or 1) with the higher contribution
# that is stored in the delta list
# genome: An array with each gene information, see Gene class
# returns the best candidate solution so far and its f value
def success_trial():
  global contribution
  D = len(contribution)
  success = False
  k = 0
  while( k<D and not success ): 
    allele, cont, trial = best_allele(k)
    if(cont != 0): success = (trial+2 >= len(contribution[k]))
    k += 1
  return success


# Gene characterization analysis trials
# genome: An array with each gene information, see Gene class
# f: Function to be optimized
# x: initial point
# evals: Maximum number of fitness evaluations
# fx: f value at point x (if provided)
def GCSA( f, evals, x, fx=None ):
  global separable, intron, contribution
  D = len(x) # Space dimension

  x, fx, evals = GCA(f, evals, x, fx) # Best solution obtained with GCA

  if(for_all(separable) or for_all(intron)): return x, fx

  failTrials = 0 if success_trial() else 1
  while(evals>=2 and failTrials<3 ):
    y, fy, evals = GCA(f, evals, bitstring(D)) # Best solution obtained with SLA
    failTrials = 0 if success_trial() else failTrials + 1
    x, y, fx, fy = pick(x,y,fx,fy)
  if(evals>0):
    y = [best_allele(k)[0] for k in range(D)] 
    fy = evaluate(f,y)
    x, y, fx, fy = pick(x, y, fx, fy)
    evals -= 1
  if(PRINT): print('Best GCA trials:', fx, x)
  return x, fx, evals


# Introns only search algorithm
# genome: An array with each gene information, see Gene class
# f: Function to be optimized
# x: initial point
# evals: Maximum number of fitness evaluations
# fx: f value in the x point, if available
def IOSA( f, evals, x, fx=None ):
  global intron
  print('Evals....',evals)
  if(not fx): fx = evaluate(f, x) # Evaluates f on x if not done
  D = len(x) # Space dimension

  introns = []
  for k in range(D):
    if(intron[k]): introns.append(k)

  if(PRINT): 
    print('Introns:', introns)
    print('Removes: (intron, remaining evaluations)')
  
  N = len(introns)
  while(evals>0 and N>0):
    j = rand.randint(0,N-1)
    k = introns[j]
    y = flip(x,k)
    fy = evaluate(f,y)
    x, y, fx, fy = pick(x, y, fx, fy)
    evals -= 1    
    # Checks if the locus is not neutral anymore and removes it
    if( C(x, y, fx, fy, k) != 0 ):
      introns.pop(j)
      N-=1
      if(PRINT): print('(', k, ',', evals, ')', sep='')
  if(PRINT): print()
  return x, fx, evals

# f: Function to be optimized
# x: initial point
# evals: Maximum number of fitness evaluations
# fx: f value at point x (if provided)
def GABO( f, evals, x, fx=None ):

  D = len(x) # Space dimension

  # Initializes component information
  init_gene(D) 

  x, fx, evals = GCSA(f, evals, x, fx) # Best solution obtained with SLA
  if(PRINT): print('******Best GCSA:', iter, fx, x)

  x, fx, evals = IOSA(f, evals, x, fx) # Best solution for 'introns'
  if(PRINT): print('******Best IOSA:', iter, fx, x) # Remove this line if printing is not required

  return x, fx