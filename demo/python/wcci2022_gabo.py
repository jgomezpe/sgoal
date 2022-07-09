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
from sgoal import evaluations
from sgoal import stats
from binary import GSC1
from binary import RMHC
from ga import basic_binary_GGA
from ga import basic_binary_SSGA
from chavela import CHAVELA
from gabo import GABO
from real import gaussian_mutation

def round(x): return (int(100*x+0.5))/100
  
# Paper reported values
def report(sgoal, fx, iter, budget, sr):
  sr *= 100/len(fx)
  avg, std = stats(fx)
  avg, std = round(avg), round(std)  
  avg_iter, std_iter = stats(iter)
  avg_iter, std_iter = round(avg_iter), round(std_iter)  
  avg_budget, std_budget = stats(budget)
  avg_budget, std_budget = round(avg_budget), round(std_budget)  
  print(sgoal, sr, avg, std, avg_iter, std_iter, avg_budget,  std_budget)
  return avg, std, avg_iter, std_iter, avg_budget, std_budget

EXP = 100 # Number of experiments

D = 120 # Bitstring length
print('=================', D, '=================')

MAX_EVALS = 100*D # Maximum number of function evaluations carried on by gabo (may require less thatn those)

#Function
testbed = [maxones, deceptive, boundedly, royalroad8, mixed] # Testbed
name = ['MaxOnes','GD3','GBD4','RR1','Mixed']
FUNCTION = 0 # Testing the deceptive function. Change the number accordingly
OPTIMUM = [D, 10*D, D, D, 47*D//20] # Optimum value of the associated test function

print('***************', name[FUNCTION], '***************')

print('SGOAL SR MAX STD ITER STD BUDGET STD')
fx = [] # Function value reached by the sgoal
iter = [] # Iter when the best value is reached by the sgoal
budget = [] # Number of fuction evaluations carried on by the sgoal
sr = 0 # Success rate of the sgoal

#GABO experiment
for i in range(EXP):
  init()
  y, fy = GABO(testbed[FUNCTION], MAX_EVALS, bitstring(D))
  fx.append(fy)
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  iter.append(best_evaluation()[1])
  budget.append(evaluations())
report('GABO', fx, iter, budget, sr)

#GSC1 experiment
fx = [] 
iter = [] 
budget = []
sr = 0
for i in range(EXP):
  init()
  y, fy = GSC1(testbed[FUNCTION], MAX_EVALS, bitstring(D))
  fx.append(fy)
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  iter.append(best_evaluation()[1])
  budget.append(MAX_EVALS)
report('GSC1', fx, iter, budget, sr)

#RMHC experiment
fx = [] 
iter = [] 
budget = []
sr = 0
for i in range(EXP):
  init()
  y, fy = RMHC(testbed[FUNCTION], MAX_EVALS, bitstring(D))
  fx.append(fy)
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  iter.append(best_evaluation()[1])
  budget.append(MAX_EVALS)
report('RMHC', fx, iter, budget, sr)

#GGA
fx = [] 
iter = [] 
budget = []
sr = 0
for i in range(EXP):
  init()
  P, fP, evals = basic_binary_SSGA(testbed[FUNCTION], MAX_EVALS, 100, D)
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  fy, k = best_evaluation()
  fx.append(fy)
  iter.append(k)
  budget.append(MAX_EVALS)
report('SSGA', fx, iter, budget, sr)

#GGA
fx = [] 
iter = [] 
budget = []
sr = 0
for i in range(EXP):
  init()
  P, fP, evals = basic_binary_GGA(testbed[FUNCTION], MAX_EVALS, 100, D)
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  fy, k = best_evaluation()
  fx.append(fy)
  iter.append(k)
  budget.append(MAX_EVALS)
report('GGA', fx, iter, budget, sr)

#CHAVELA
fx = [] 
iter = [] 
budget = []
sr = 0
for i in range(EXP):
  init()
  P, fP, evals, rates = CHAVELA(testbed[FUNCTION], MAX_EVALS, [simple_crossover, bit_mutation, transposition], bitstring_population(100,D))
  sr += 1 if success_evaluation(OPTIMUM[FUNCTION]) != -1 else 0
  fy, k = best_evaluation()
  fx.append(fy)
  iter.append(k)
  budget.append(MAX_EVALS)
report('CHAVELA', fx, iter, budget, sr)

