# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Classical Genetic Algorithms
from sgoal.core import evaluate_population
from sgoal.core import evaluate
from sgoal.core import best
from sgoal.core import tracing
from sgoal.core import tournament
from sgoal.util import randbool
from sgoal.core import simple_crossover
from sgoal.binary import bit_mutation_probability
from sgoal.binary import bitstring_population

# Tracing the population's fitness evolution
FP = []

# Gets the genetic algorithm's traced information
def ga_trace():
  global FP
  return FP
  
############### Generational Genetic Algorithm - GGA ################
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# selection: Parents selection mechanism
# xr: Crossover rate
# xover: Crossover operator
# mutation: Mutation operator
# P: Initial population
# fP: Fitness value for each individual of the population, if provided.
def GGA(f, evals, selection, xr, xover, mutation, P, fP=None):
  global FP # Tracing information
  n = len(P)
  if( evals>=n and not fP ):
    fP = evaluate_population(f, P)
    evals -= n
  if( tracing() ): FP = [fP]
  while( evals>=2 ):    
    Q = []
    fQ = []
    for i in range(0,n,2):
      idx1, idx2 = selection(fP, 2)
      if evals>=2 and randbool(xr):
        a, b = xover(P[idx1], P[idx2])
        a = mutation(a)
        b = mutation(b)
        Q.append(a)
        Q.append(b)
        fQ.append(evaluate(f,a))
        fQ.append(evaluate(f,b))
        evals -= 2
      else:
        Q.append( P[idx1] )
        Q.append( P[idx2] )
        fQ.append(fP[idx1] )
        fQ.append(fP[idx2] )
    P = Q
    fP = fQ
    if( tracing() ): FP.append(fP)
  return P, fP, evals

############### Steady State Genetic Algorithm - GGA ################
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# selection: Parents selection mechanism
# xover: Crossover operator
# mutation: Mutation operator
# P: Initial population
# fP: Fitness value for each individual of the population, if provided.
def SSGA(f, evals, selection, xover, mutation, P, fP=None):
  global FP # Tracing information
  n = len(P)
  if( evals>=n and not fP ):
    fP = evaluate_population(f, P)
    evals -= n
  if( tracing() ): FP = [fP]
  i = 0
  while( evals>=2 ):  
    idx1, idx2 = selection(fP, 2)
    p1, p2 = P[idx1], P[idx2]
    a, b = xover(p1, p2)
    sol = [p1, p2, mutation(a), mutation(b)]
    fsol = [fP[idx1], fP[idx2], evaluate(f,sol[2]), evaluate(f,sol[3])]
    evals -= 2
    fP[idx1], k = best(fsol)
    P[idx1] = sol[k]
    sol.pop(k)
    fsol.pop(k)
    fP[idx2], k = best(fsol)
    P[idx2] = sol[k]
    i += 2
    if(i==n and tracing()):
      FP.append(fP.copy())
      i=0
  return P, fP, evals

############### Basic Binary Generational Genetic Algorithm ################
# Uses simple crossover, bit_mutation, and tournament mechanism as parents selection
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# Mu: Population's size
# L: Bitstring's length
# xr: Crossover rate (0.7 by default)
# mr: Mutation rate (1/L by default)
def basic_binary_GGA(f, evals, Mu, L, xr=0.7, mr=None):
  mr = 1.0/L if not mr else mr
  mutation = lambda x: bit_mutation_probability(x,mr)
  P = bitstring_population(Mu, L)
  return GGA(f, evals, tournament, xr, simple_crossover, mutation, P)

############### Basic Binary Steady State Genetic Algorithm ################
# Uses simple crossover, bit_mutation, and tournament mechanism as parents selection
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# Mu: Population's size
# L: Bitstring's length
# mr: Mutation rate (1/L by default)
def basic_binary_SSGA(f, evals, Mu, L, mr=None):
  mr = 1.0/L if not mr else mr
  mutation = lambda x: bit_mutation_probability(x,mr)
  P = bitstring_population(Mu, L)
  return SSGA(f, evals, tournament, simple_crossover, mutation, P)