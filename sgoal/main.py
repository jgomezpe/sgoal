from binary import BitArrayProblem
from binary import BitArrayCHAVELA
from binary import BitArrayR1_5
from binary import BitArrayGGA
from binary import BitArraySSGA
from binary import GS1
from binary import GSC1
from binary import RMHC
from binary import BitArrayHC
from gabo import GABO
from gabo2 import GABO2
from sgoal import experiment

#Main program
R = 10
D = 120
EVALS = 100*D
problem = BitArrayProblem('mixed', D)
sgoal = [GABO2, GABO, BitArrayHC, RMHC, BitArrayR1_5, GS1, GSC1, BitArrayGGA, BitArraySSGA, BitArrayCHAVELA]

for s in sgoal:
  f, evals, sr = experiment(s, problem, EVALS, R)
  print( sum(f)/R, sum(evals)/R, sr) 

