import pandas as pd
import random as rand
from sgoal.core import randbool
from sgoal.core import rec
from sgoal.binary import Binary
from sgoal.binary import flip
from sgoal.binary import multiflip
from sgoal.core import PROBLEM
from sgoal.core import SPSGoal
from sgoal.core import simplereplace
from sgoal.core import next

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

def maxcut(x, W):
  s = 0
  for w in W:
    if(x[w[0]]+x[w[1]]==1):
      s += w[2]
  return s

def fastmaxcut(y, k, x, fx, sgoal):
  REL, WREL, M = sgoal['REL'], sgoal['WREL'], len(sgoal['W'])
  count = 0
  n = len(x)
  if(len(k)==n or len(k)==0): return fx
  f = fx
  for i in k:
    for j in range(len(REL[i])): 
      r = REL[i][j]
      if(x[i]==x[r] and y[i]!=y[r]):
        count += 1
        f += WREL[i][j]
      elif(x[i]!=x[r] and y[i]==y[r]):
        count += 1
        f -= WREL[i][j]
  sgoal['delta'] = 3*count/M
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

def bitmutation(x, fx, sgoal):
  p = 1/len(x)
  k = []
  for i in range(len(x)):
    if(randbool(p)):
      k.append(i)
  if(len(k)>0): y = multiflip(x, k)
  else: y = x.copy()
  fy = fastmaxcut(y, k, x, fx, sgoal)
  return y, fy

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

def readG(T):
  return read(Gdimension(T), 'https://web.stanford.edu/~yyye/yyye/Gset/G' + str(T))


def GProblem(REL, WREL, W, EVALS):
  D = len(REL)
  space = Binary(D)
  problem = PROBLEM('max', lambda x: maxcut(x,W), space, EVALS)
  problem['REL'] = REL
  problem['WREL'] = WREL
  problem['W'] = W
  problem['flip'] = lambda x, fx, k: sflip(x, fx, k, problem)
  problem['multiflip'] = lambda x, fx, k: mflip(x, fx, k, problem)
  return problem

def maxcutHC(problem):
  if('variation' not in problem): problem['variation'] = lambda x, fx : singlebitmutation(x, fx, problem)
  if('replace' not in problem): problem['replace'] = lambda x,fx, y, fy : simplereplace(x, fx, y, fy, problem)
  problem['next'] = lambda x, fx : next(x, fx, problem)
  return SPSGoal(problem)