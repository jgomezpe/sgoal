from binary import maxones
from binary import deceptive
from binary import boundedly
from binary import royalroad8
from binary import mixed
from binary import bitstring
from binary import bit_mutation
from binary import bitstring_population
from sgoal import simple_crossover
from sgoal import transposition
from sgoal import init
from sgoal import best_evaluation
from sgoal import success_evaluation
from sgoal import stats
from sgoal import MAXIMIZE
from binary import GSC1
from binary import RMHC
from gabo import GABO

# Paper reported values
def report(sgoal, fx, iter, sr):
  sr *= 100/len(fx)
  avg, std = stats(fx)
  avg_iter, std_iter = stats(iter)
  print("Results for", sgoal)
  print("success rate: ", sr, "%", sep='')
  print("best: ", avg, "+/-", std, sep='')
  print("iter: ", avg_iter, "+/-", std_iter, sep='')
  print('*********************************')
  return avg, std, avg_iter, std_iter

# Global variables
testbed = [maxones, deceptive, boundedly, royalroad8, mixed] # Testbed
name = ['MaxOnes','GD3','GBD4','RR1','Mixed']
LENGTH = [100, 30, 40, 64, 100] # Bitstring length
MAX_EVALS = [205, 500, 1000, 6400, 5000] # Maximum number of function evaluations
OPTIMUM = [LENGTH[0], 10*LENGTH[1], LENGTH[2], LENGTH[3], 47*LENGTH[4]//20] # Optimum value of the associated test function
EXP = 100 # Number of experiments

for FUNCTION in range(len(testbed)): # Testing the deceptive function. Change the number accordingly
  print('***************', name[FUNCTION], '***************')
  fx = [] # Function value reached by the sgoal
  iter = [] # Iter when the best value is reached by the sgoal
  sr = 0 # Success rate of the sgoal

  #GABO experiment
  for i in range(EXP):
    init()
    y, fy = GABO(testbed[FUNCTION], MAX_EVALS[FUNCTION], bitstring(LENGTH[FUNCTION]))
    fx.append(fy)
    sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
    iter.append(best_evaluation()[1])
  report('GABO', fx, iter, sr)

  #GSC1 experiment
  fx = [] 
  iter = [] 
  sr = 0
  for i in range(EXP):
    init()
    y, fy = GSC1(testbed[FUNCTION], MAX_EVALS[FUNCTION], bitstring(LENGTH[FUNCTION]))
    fx.append(fy)
    sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
    iter.append(best_evaluation()[1])
  report('GSC1', fx, iter, sr)

  #RMHC experiment
  fx = [] 
  iter = [] 
  sr = 0
  for i in range(EXP):
    init()
    y, fy = RMHC(testbed[FUNCTION], MAX_EVALS[FUNCTION], bitstring(LENGTH[FUNCTION]))
    fx.append(fy)
    sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
    iter.append(best_evaluation()[1])
  report('RMHC', fx, iter, sr)
