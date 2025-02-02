from sgoal.binary import BitArrayProblem
from sgoal.binary import BitArrayCHAVELA
from sgoal.binary import BitArrayR1_5
from sgoal.binary import BitArrayGGA
from sgoal.binary import BitArraySSGA
from sgoal.binary import GS1
from sgoal.binary import GSC1
from sgoal.binary import RMHC
from sgoal.binary import BitArrayHC
from sgoal.gabo import GABO
from sgoal.gabo2 import GABO2
from sgoal.core import experiment

#Main program
R = 10
D = 120
EVALS = 100*D
problem = BitArrayProblem('Mixed', D)
core = [GABO2, GABO, BitArrayHC, RMHC, BitArrayR1_5, GS1, GSC1, BitArrayGGA, BitArraySSGA, BitArrayCHAVELA]

for s in core:
  f, evals, sr = experiment(s, problem, EVALS, R)
  print( sum(f)/R, sum(evals)/R, sr) 

