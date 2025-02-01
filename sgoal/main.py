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
from gabosfm import GABOFSI


#Main program
D = 120
problem = binaryproblem('mixed', D)

s = []
for i in range(30):
    sgoal = GABOFSI(problem)
    r = sgoal.run(100*D, True)
    print(r['groups'])
    s.append(r['f'])
#print(r['trace'])
#print(r['besttrace'])
    print(i, r['evals'], r['f'])
#print(r['poptrace'])
print(sum(s)/len(s))

s = []
for i in range(30):
    sgoal = GABO(problem)
    r = sgoal.run(100*D, True)
    s.append(r['f'])
#print(r['trace'])
#print(r['besttrace'])
    print(i, r['evals'], r['f'])
#print(r['poptrace'])
print(sum(s)/len(s))

def dummu():
    for i in range(1):
        sgoal = bitCHAVELA(problem)
        r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
        print(i, r['evals'], r['f'])
    #print(r['poptrace'])

    sgoal = bitR1_5(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])

    sgoal = HC(problem, bitmutation)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    sgoal = RMHC(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    sgoal = DGS1(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    sgoal = DGSC1(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    sgoal = bitGGA(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    #print(r['poptrace'])
    sgoal = bitSSGA(problem)
    r = sgoal.run(100*D, True)
    #print(r['trace'])
    #print(r['besttrace'])
    print(r['evals'], r['f'])
    #print(r['poptrace'])
