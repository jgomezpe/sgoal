import math
import random as rand
from sgoal.core import randbool
from sgoal.core import Space

# Fixed Length RealSpace
# min: Min values
# max: Max values
class RealSpace(Space):
  def __init__(self, min, max):
    self.min = min
    self.max = max
    if(isinstance(min , list)):
       self.length = [max[i]-min[i] for i in range(len(min))]
    else:
      self.length = max - min

  # Produces one random BitArray of length D
  def getone(self):
    if(isinstance(self.min , list)):
      return [self.min[i] + rand.random()*self.length[i] for i in range(len(self.min))]
    return self.min + rand.random()*self.length

  # Determines if a candidate solution is feasible
  def feasible(self, x):
    if(isinstance(self.min , list)):
      for i in range(len(x)):
        if(not (self.min[i] <= x[i] <= self.max[i])):
            return False
    return self.min <= x <= self.max


# Guassian number generation
def gaussian(sigma=1.0): return rand.gauss(sigma)

# Power law number generation
def powerlaw(a=-2.0):
  if a==-2.0: return 1/(1-rand.random())
  else: return (1-rand.random())**(1/(a+1))

############### VARIATION OPERATIONS ################
# Applies a mutation to each real component with probability p
def intensity_mutation(x, mutation, p):
  y = x.copy()
  for i in range(len(y)):
    if( randbool(p) ): y[i] = mutation(y[i]) 
  return y

# Gaussian Mutation
def gaussianMutation( x, sigma, space ):
  if( isinstance(x, list) ):
    y = x.copy()
    i = rand.randint(0,len(x)-1)
    y[i] += gaussian(sigma) # = [z + gaussian(sigma) for z in x]
    return y if space.feasible(y) else x
  y = x + gaussian(sigma)
  return y if space.feasible(y) else x

# A lambda version for using in SGoals requiring one argument variations
def lambdaGaussianMutation(sigma, space):
   return lambda x: gaussianMutation(x, sigma, space)

# A lambda version for using in SGoals requiring one argument variations
def lambdaOneGaussianMutation(space):
   return lambda x, sigma: gaussianMutation(x, sigma, space)

# Powerlaw Mutation
def powerlawMutation( x, space ):
  if( isinstance(x, list) ):
    y = [z + (1 if randbool() else -1)*powerlaw() for z in x]
    return y if space.feasible(y) else x
  y = x + powerlaw()
  return y if space.feasible(y) else x

# A lambda version for using in SGoals requiring one argument variations
def lambdaPowerlawMutation(space):
   return lambda x: powerlawMutation(x, space)
   
#################### TEST FUNCTIONS ##############
# Sphere function
def sphere_1(x):
   return x*x

def sphere(x):
   s = 0.0
   for y in x:
      s += y*y
   return s

# Rastrigin function as proposed by Rastrigin, L. A. in "Systems of extremal control." Mir, Moscow (1974).
def rastrigin_1( x ):
	return x*x - 10.0*math.cos(2.0*math.pi*x)

def rastrigin( x ):
  f = 0.0
  for c in x:
    f += rastrigin_1(c)
  return 10.0*len(x) + f

# Schwefel Function
def schwefel_1( x ):
	return -x * math.sin(math.sqrt(math.abs(x)))

def schwefel( x ):
  f = 0.0
  for c in x:
    f += schwefel(c)
  return (418.9829101*len(x) + f)

# Griewangk function
def griewank( x ):
  sum = 0.0
  prod = 1.0
  for i in range(len(x)):
    sum += x[i]*x[i]/4000.0
    prod *= math.cos(x[i]/math.sqrt(i+1.0))
  return (1.0 + sum - prod)

# Rosenbrock Saddle Function
def rosenbrock_saddle_2( x1, x2 ):
	y = x1*x1 - x2
	return (100.0*y*y + (1.0-x1)*(1.0-x1))

def rosenbrock_saddle( x ):
  f = 0.0
  for i in range(len(x)-1):
    f += rosenbrock_saddle_2( x[i], x[i+1] )
  return f

##################### TEST PROBLEMS ####################
def RealProblem(f, D):
  if(f=='Rastrigin'):
    return {'f':rastrigin, 'space': RealSpace([-5.12 for i in range(D)],[5.12 for i in range(D)]), 'optimum':0.0, 'type':'min'}
  if(f=='Schwefel'):
    return {'f':schwefel, 'space': RealSpace([-500.0 for i in range(D)],[500.0 for i in range(D)]), 'optimum':0.0, 'type':'min'}
  if(f=='Griewank'):
    return {'f':griewank, 'space': RealSpace([-600.0 for i in range(D)],[600.0 for i in range(D)]), 'optimum':0.0, 'type':'min'}
  if(f=='Rosenbrock'):
    return {'f':rosenbrock_saddle, 'space': RealSpace([-2.048 for i in range(D)],[2.048 for i in range(D)]), 'optimum':0.0, 'type':'min'}
  if(f=='Sphere'):
    return {'f':sphere, 'space':  RealSpace([-5.12 for i in range(D)],[5.12 for i in range(D)]), 'optimum':0.0, 'type':'min'}
  return {'f':sphere, 'space':  RealSpace([-5.12 for i in range(D)],[5.12 for i in range(D)]), 'optimum':0.0, 'type':'min'}