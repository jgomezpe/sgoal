# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Classical Genetic Algorithms
from sgoal import SGoal
from sgoal import randbool
from sgoal import best

############### Generational Genetic Algorithm - GGA ################
# problem: Problem to be solved
# config: GGA parameters
class GGA(SGoal):
  def __init__(self, problem, config):
    SGoal.__init__(self, problem)
    self.selection = config['selection']
    self.mutation = config['mutation']
    self.xr = config['xr']
    self.xover = config['xover']
    self.N = config['N']
    if(self.N%2==1): 
      self.N+=1
    self.poptrace = []

  def next(self, P, fP):
    Q = []
    fQ = []
    for i in range(0,self.N,2):
      idx1, idx2 = self.selection(fP, 2, self.minimize)
      if self.caneval(2) and randbool(self.xr):
        a, b = self.xover(P[idx1], P[idx2])
        a = self.mutation(a)
        b = self.mutation(b)
        Q.append(a)
        Q.append(b)
        fQ.append(self.evalone(a))
        if(not self.stop()): 
          fQ.append(self.evalone(b) )
      else:
        Q.append( P[idx1] )
        Q.append( P[idx2] )
        fQ.append(fP[idx1] )
        fQ.append(fP[idx2] )
    P = Q
    fP = fQ
    return P, fP

############### Steady State Genetic Algorithm - GGA ################
class SSGA(GGA):
  def __init__(self, problem, config):
    GGA.__init__(self, problem, config)

  def next(self, P, fP):
    i=0
    while(i<self.N and self.caneval(2)):
      idx1, idx2 = self.selection(fP, 2, self.minimize)
      p1, p2 = P[idx1], P[idx2]
      a, b = self.xover(p1, p2)
      a, b = self.mutation(a), self.mutation(b)
      fa, fb = self.evalone(a), self.evalone(b)
      k = best(fP, not self.minimize)
      P[k] = a
      fP[k] = fa
      k = best(fP, not self.minimize)
      P[k] = b
      fP[k] = fb
      i += 2
    return P, fP
