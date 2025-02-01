# Copyright (c)
# Author: Jonatan Gomez 
# E-mail: jgomezpe@unal.edu.co
# All rights reserved.
# Classical Hill Climbing Algorithm with neutral mutations
from sgoal.core import SGoal

# Classical Hill Climbing Algorithm with neutral mutations
class HC(SGoal):
  def __init__(self, problem, variation):
    SGoal.__init__(self, problem)
    self.variation = variation

  def next(self, P, fP):
    y = self.variation(P)
    fy = self.eval(y)
    P, fP, y, fy = self.pick(P, fP, y, fy)
    return P, fP