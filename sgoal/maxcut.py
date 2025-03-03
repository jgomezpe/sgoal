import pandas as pd
import random as rand
from sgoal.core import randbool
from sgoal.core import rec
from sgoal.binary import Binary
from sgoal.binary import flip
from sgoal.binary import multiflip
from sgoal.core import PROBLEM
from sgoal.core import VRSGoal
from sgoal.core import simplereplace
from sgoal.core import next
from sgoal.core import variation11
from sgoal.es import Rule1_5_T
from sgoal.ga import SSGA_T
from sgoal.ga import GGA_T
from sgoal.chavela import apply
from sgoal.chavela import CHAVELA_T
from sgoal.chavela import CHAVELA1_T

# Read a maxcut problem matrix
def read(D, url):
  REL = [[] for i in range(D)]
  WREL = [[] for i in range(D)]
  W = []
  data = pd.read_csv(url, sep=' ', skiprows=1, header=None, names=None)
  for w in data.values:
    a = w[0]-1
    b = w[1]-1 
    REL[a].append(b)
    REL[b].append(a)
    WREL[a].append(w[2])
    WREL[b].append(w[2])
    W.append([a,b,w[2]])
  return REL, WREL, W

########## Standard version of maxcut problem ##########
# Standard function evaluation
def maxcut(x, W):
  s = 0
  for w in W:
    if(x[w[0]] != x[w[1]]):
      s += w[2]
  return s

# Standard Maxcut Problem
def MaxCutProblem(REL, WREL, W, EVALS):
  D = len(REL)
  space = Binary(D)
  problem = PROBLEM('max', lambda x: maxcut(x,W), space, EVALS)
  problem['REL'] = REL
  problem['WREL'] = WREL
  problem['W'] = W
  return problem


########## Fast version of maxcut problem ##########
def fastmaxcut(y, k, x, fx, sgoal):
  REL, WREL, M = sgoal['REL'], sgoal['WREL'], len(sgoal['W'])
  count = 0
  n = len(x)
  if(len(k)==n or len(k)==0): return fx
  delta = 0
  for i in k:
    if(x[i]==y[i]): print('what')
    for j in range(len(REL[i])): 
      r = REL[i][j]
      if(x[r] == y[r]):
        if(x[i]==x[r]):
          count += 1
          delta += WREL[i][j]
        else:
          count += 1
          delta -= WREL[i][j]
  f = fx + delta
  #sgoal['delta'] = 3*count/M
  rec(y, f, sgoal)
  return f 

def sflip(x, fx, k, sgoal):
  y = flip(x, k)
  fy = fastmaxcut(y, [k], x, fx, sgoal)
  return y, fy

def mflip(x, fx, k, sgoal):
  y = multiflip(x, k)
  fy = fastmaxcut(y, k, x, fx, sgoal)
  return y, fy

def singlebitmutation(x, fx, sgoal):
  k = rand.randint(0,len(x)-1)
  y = flip(x, k)
  fy = fastmaxcut(y, [k], x, fx, sgoal)
  return y, fy

def bitmutationprob(x, fx, p, sgoal):
  k = []
  for i in range(len(x)):
    if(randbool(p)):
      k.append(i)
  if(len(k)>0): y = multiflip(x, k)
  else: y = x.copy()
  fy = fastmaxcut(y, k, x, fx, sgoal)
  return y, fy

def bitmutation(x, fx, sgoal):
  return bitmutationprob(x,fx,1/len(x),sgoal)

# Simple crossover.
def xover( x1, x2 ):
  n = len(x1)
  p = rand.randint(1,n-1)
  y1 = x1[0:p] + x2[p:n]
  y2 = x2[0:p] + x1[p:n]
  r1 = range(p,n)
  r2 = range(0,p)
  if(randbool()):
    y1, y2 = y2, y1
    r1, r2 = r2, r1
  return y1, y2, r1, r2

def simplexover( x1, fx1, x2, fx2, sgoal ):
  y1, y2, r1, r2 = xover(x1, x2)
  k1 = []
  for i in r1:
    if(x1[i] != y1[i]):
      k1.append(i)
  k2 = []
  for i in r2:
    if(x2[i] != y2[i]):
      k2.append(i)
  fy1 = fastmaxcut(y1, k1, x1, fx1, sgoal)
  fy2 = fastmaxcut(y2, k2, x2, fx2, sgoal)
  return y1, fy1, y2, fy2

def simplexover1( x, fx, sgoal ):
  P, fP, selection = sgoal['P'], sgoal['fP'], sgoal['selection']
  idx = selection(fP, 1)[0]
  x2 = P[idx]
  y, y2, r, r2 = xover(x, x2)
  k = []
  for i in r:
    if(x[i] != y[i]):
      k.append(i)
  fy = fastmaxcut(y, k, x, fx, sgoal)
  return y, fy


# Transposition
def transposition( x, fx, sgoal ):
  y = x.copy()
  start = rand.randint(0,len(x)-1)
  end = rand.randint(0,len(x)-1)
  if start>end: start, end = end, start
  r = range(start, end)
  while start<end:
    y[start], y[end] = y[end], y[start]
    start += 1
    end -= 1
  k = []
  for i in r:
    if(x[i] != y[i]):
      k.append(i)
  fy = fastmaxcut(y, k, x, fx, sgoal)
  return y, fy


##################### SGOALs ###########################
# Classical Hill Climbing Algorithm with neutral mutation for BitArray problems. Uses bitmutation as variation operator
# problem: Problem to solve
def HC(problem): 
  if( 'variation' not in problem ): problem['variation'] = bitmutation 
  return VRSGoal(problem)

# The HC algorithm suggested by Richard Palmer, that Forrest and
# Mitchell named as "random mutation hill-climbing" (RMHC), see
# M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?”
# Santa Fe Institute, Working Papers, 01 1993.
# problem: Problem to solve
def RMHC(problem): 
  problem['variation'] = lambda x, fx : singlebitmutation(x, fx, problem)
  return VRSGoal(problem)

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, for BitArray
# problem: Problem to solve
def setprob(problem):
  problem['variation'] = lambda x, fx: variation11(x, fx, lambda y: bitmutationprob(y, problem['parameter']), problem)
  
# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
def Rule1_5(problem):
  D = problem['D']
  if( 'parameter' not in problem ): problem['parameter'] = 1/D
  if( 'variation' not in problem ): 
    problem['setparameter'] = lambda : setprob(problem)
    setprob(problem)
  if( 'G' not in problem ): problem['G'] = D
  return Rule1_5_T(problem)

############### Generational Genetic Algorithm - SSGA ################
def bmutation(sgoal):
  if('mutation' not in sgoal): sgoal['mutation'] = bitmutation
  return sgoal

def GGA(problem):
  return GGA_T(bmutation(problem))

############### Steady State Genetic Algorithm - SSGA ################
# problem: Problem to solve
def SSGA(problem):
  return SSGA_T(bmutation(problem))

# Standard CHAVELA for Binary problems. Uses bitmutation, simplexover, and transposition as operators
def CHAVELA(problem):
  mutation = lambda x, fx: bitmutation(x, fx, problem)
  xover = lambda x, fx: simplexover1(x, fx, problem)
  transp = lambda x, fx: transposition( x, fx, problem)
  if( 'operators' not in problem ): problem['operators'] = [mutation, transp, xover]
  return CHAVELA_T(problem)

# Standard CHAVELA1 for Binary problems
#def CHAVELA1(problem):
#  if( 'operators' not in problem ): problem['operators'] = [powerlawmutation, singlebitmutation]
#  return CHAVELA1_T(problem) 



def FastMaxCutProblem(REL, WREL, W, EVALS):
  problem = MaxCutProblem(REL, WREL, W, EVALS)
  problem['flip'] = lambda x, fx, k: sflip(x, fx, k, problem)
  problem['multiflip'] = lambda x, fx, k: mflip(x, fx, k, problem)
  return problem

########### Test bed #############
### G data set from https://web.stanford.edu/~yyye/yyye/Gset/ 
# Read a maxcut problem matrix from the G data set 
def readG(T):
  return read(Gdimension(T), 'https://web.stanford.edu/~yyye/yyye/Gset/G' + str(T))

def Gdimension(T):
  D=800
  if( T<=21):
    D=800
  elif( 21<T<=42 ):
    D=2000
  elif(T<=47):
    D=1000
  elif(T<=59):
    D=5000
  elif(T<=64):
    D=7000
  elif(T==65):
    D=8000
  elif(T==66):
    D=9000
  elif(T<=72):
    D=10000
  elif(T==77):
    D=14000
  else:
    D=20000
  return D
