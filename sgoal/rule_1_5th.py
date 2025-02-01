# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
from sgoal import SGoal

# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# G: Checking evaluations for adapting the mutation rate
# a: Scaling factor of the mutation rate
# mutation: A mutation operation with the possibility of setting its rate
# mr: Initial mutation rate
# x: Initial point
# fx: f value at point x (if provided)
class Rule_1_5th(SGoal):
  def __init__(self, problem, config):
    SGoal.__init__(self, problem)
    self.G = config['G']
    self.a = config['a']
    self.mr = config['mr']
    self.mutation = config['mutation']
    self.I = 0
    self.Gs = 0
    
  def init(self, MAXEVALS, TRACE):
    self.I = 1 
    self.Gs = 1
    return SGoal.init(self, MAXEVALS, TRACE)

  def next(self, P, fP):
    y = self.mutation(P, self.mr)
    fy = self.evalone(y)
    w = P
    P, fP, y, fy = self.pick( P, fP, y, fy )
    self.I += 1
    if( w==y ): self.Gs += 1
    if( self.I==self.G ):
      p = self.Gs/self.G
      if( p > 0.2 ): self.mr /= self.a
      elif( p < 0.2 ): self.mr *= self.a
      self.Gs = 0
      self.I = 0    
    return P, fP