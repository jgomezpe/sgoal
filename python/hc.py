# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Classical Hill Climbing Algorithm with neutral mutations
from sgoal import pick
from sgoal import evaluate

# Classical Hill Climbing Algorithm with neutral mutations
# f: Function to be optimized
# evals: Maximum number of fitness evaluations
# variation: Variation operation
# x: Initial point
# fx: f value at point x (if provided)
def HC( f, evals, variation, x, fx=None):
  if(evals>0 and not fx): 
    fx = evaluate(f, x)
    evals-=1
  for i in range(evals):
    y = variation(x)
    fy = evaluate(f,y)
    x, y, fx, fy = pick( x, y, fx, fy )
  return x, fx