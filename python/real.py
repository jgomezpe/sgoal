import math
import numpy
import random as rand
from sgoal import randbool 

def rastrigin_1( x ):
	return x*x - 10.0*math.cos(2.0*math.pi*x)

def rastrigin( x ):
  f = 0.0
  for c in x:
    f += rastrigin_1(c)
  return 10.0*len(x) + f

def schwefel_1( x ):
	return -x * math.sin(math.sqrt(math.abs(x)))

def schwefel( x ):
  f = 0.0
  for c in x:
    f += schwefel(c)
  return (418.9829101*len(x) + f)

def griewangk( x ):
  sum = 0.0
  prod = 1.0
  for i in range(len(x)):
    sum += x[i]*x[i]/4000.0
    prod *= math.cos(x[i]/math.sqrt(i+1.0))
  return (1.0 + sum - prod)

def rosenbrock_saddle_2( x1, x2 ):
	y = x1*x1 - x2
	return (100.0*y*y + (1.0-x1)*(1.0-x1))

def rosenbrock_saddle( x ):
  f = 0.0
  for i in range(len(x)-1):
    f += rosenbrock_saddle_2( x[i], x[i+1] )
  return f

def gaussian(sigma=1.0): return numpy.random.normal(scale=sigma)[0]
  
def powerlaw(a=-2.0):
  if a==-2.0: return 1/(1-rand.random())
  else: return (1-rand.random())**(1/(a+1))

# Applies a mutation to each real component with probability p
def intensity_mutation(x, mutation, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): y[i] = mutation(y[i]) 
  return y

def gaussian_mutation( x, sigma ):
  d = numpy.random.normal(scale=sigma,size=len(x))
  print(d)
  y = x.copy()
  for i in range(len(x)): y[i] += d[i]
  return y

