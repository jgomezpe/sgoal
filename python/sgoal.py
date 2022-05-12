import random as rand
from inspect import signature

############### USEFUL FUNCTIONS ################
# Generates a boolean value according to probability p ( True with probability p, False otherwise )
def randbool(p=0.5):
  return (rand.random() < p)

# A permutation of n elements
def permutation(n):
  x = [i for i in range(0,n)]
  rand.shuffle(x)
  return x

# Normalizes a vector of weights
def normalize(weight):
  c = 0
  for w in weight: c += w
  return [ w/c for w in weight ]

####################### Global variables ########################
MAXIMIZE = True # Set to false if minimizing is required
TRACE = True # Tracing f value evolution
PRINT = True # Printing intermediate information

candidate = [] # Candidate solutions 
generated = [] # f value of each generated candidate solution

# Initializes global variables. 
# WARNING: Before running any sgoal, call this method if you are tracing the method
def init():
  global generated, candidate
  generated = []
  candidate = []

# Evaluates the function on a given candidate solution and stores information in the global generated list if tracing
def evaluate(f, x):
  global generated, TRACE
  fx = f(x)
  if( TRACE ): generated.append(fx)
  return fx

# Evaluates the function on a given population of candidate solutions and stores information in the global generated list if tracing
def evaluate_population(f, P):
  return [evaluate(f,x) for x in P]

# Sorts two candidate solutions according to their f values
# If x is the current solution and y is a new candidate solution obtained by a 
# mutation operation, this method picks y if it has equal or higher f value 
# (neutral mutations are allowed) otherwise returns x (hill climbing replacement)
def pick( x, y, fx, fy ):
  global MAXIMIZE
  if((MAXIMIZE and fy>=fx) or (not MAXIMIZE and fy<=fx)): 
    return y, x, fy, fx
  return x, y, fx, fy

# Returns the position of a value. Returns -1 if not found.
def index(list, value):
  for i in range(list):
    if(list[i]==value): return i
  return -1

# Returns the number of fitness's evaluation required by the sgoal for generating the optimal solution. Returns -1 if not found.
def success_evaluation(optimum):
  global generated
  return index(generated, optimum)

# Returns the best and the index of the best
def best(results):
  k = 0
  for i in range(1, len(results)):
    if((MAXIMIZE and results[i]>results[k]) or (not MAXIMIZE and results[i]<results[k])): 
      k = i
  return results[k], k

# Returns the number of fitness's evaluation required by the sgoal for generating the best solution.
def best_evaluation():
  global generated
  return best(generated)

############### VARIATION OPERATIONS ################
# Simple crossover.
def simple_crossover( x1, x2 ):
  p = rand.randint(1,len(x1)-1)
  y1 = x1[0:p] + x2[p:len(x2)]
  y2 = x2[0:p] + x1[p:len(x1)]
  return y1, y2

# Transposition
def transposition(x):
  x = x.copy()
  start = rand.randint(0,len(x)-1)
  end = rand.randint(0,len(x)-1)
  if start>end: start, end = end, start
  while start<end:
    x[start], x[end] = x[end], x[start]
    start += 1
    end -= 1
  return x

############### SELECTION MECHANISMS ################
# Uniform selection. Pick n elements (indices) with the same probability
def uniform(M, n):
  M -= 1
  return [rand.randint(0,M) for i in range(0,n)]


# Tournament selection. Selects n individuals. 
# For each individual to be selected this method picks 4 individuals and returns the best of those
def tournament1(quality):
  global MAXIMIZE
  m = 4 #Tournament's size
  candidate = uniform(len(quality),m)
  x = 0
  for k in range(1,m):
    if( (MAXIMIZE and quality[candidate[k]]>quality[candidate[x]]) or (not MAXIMIZE and quality[candidate[k]]<quality[candidate[x]])):
      x = k
  return candidate[x]

def tournament(quality, n):
  return [ tournament1(quality) for i in range(0,n) ]

# Weighted selection: p_i is the probability of selecting element i
def weighted(p):
  y = rand.random()
  k=0
  while k<len(p) and y>=p[k]:
    y -= p[k]
    k+=1
  return k

# Arity of a function f
def arity(f):
  return len(signature(f).parameters)
