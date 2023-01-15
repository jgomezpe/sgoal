import random as rand

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

# Statistical information 
def stats(list):
  avg = 0.0
  for x in list: avg += x
  avg /= len(list)
  std = 0.0
  for x in list: std += (x-avg)**2
  std = (std/len(list))**0.5
  return avg, std
