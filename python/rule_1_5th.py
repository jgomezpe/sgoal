# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see
# Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. 
# Natural Computing. 1. 3-52. 10.1023/A:1015059928466. 
from sgoal import pick
from sgoal import evaluate

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
def rule_1_5th( f, evals, G, a, mutation, mr, x, fx=None):
  I=0
  if(evals>0 and not fx):
    fx = evaluate(f, x)
    evals-=1
    I = 1
  Gs = 1
  for i in range(evals):
    y = mutation(x, mr)
    fy = evaluate(f,y)
    w = x
    x, y, fx, fy = pick( x, y, fx, fy )
    I += 1
    if( w==y ): Gs += 1
    if( I==G ):
      P = Gs/G
      if( P > 0.2 ): mr /= a
      elif( P < 0.2 ): mr *= a
      Gs = 0
      I = 0    
  return x, fx