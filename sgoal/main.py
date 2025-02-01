from binary import binaryproblem
from binary import bitCHAVELA
from binary import bitR1_5
from binary import bitGGA
from binary import bitSSGA
from binary import DGS1
from binary import DGSC1
from binary import RMHC
from binary import bitmutation
from hc import HC
from gabo import GABO


#Main program
D = 120
problem = binaryproblem('mixed', D)

sgoal = GABO(problem)
r = sgoal.run(100*D, True)
print(i, r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])
#print(r['poptrace'])

sgoal = bitCHAVELA(problem)
r = sgoal.run(100*D, True)
print(i, r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])
#print(r['poptrace'])

sgoal = bitR1_5(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])

sgoal = HC(problem, bitmutation)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])

sgoal = RMHC(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])

sgoal = DGS1(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])

sgoal = DGSC1(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])

sgoal = bitGGA(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])
#print(r['poptrace'])

sgoal = bitSSGA(problem)
r = sgoal.run(100*D, True)
print(r['evals'], r['f'])
#print(r['trace'])
#print(r['besttrace'])
#print(r['poptrace'])
