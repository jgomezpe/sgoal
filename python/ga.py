from sgoal import evaluate_population
from sgoal import evaluate
from sgoal import best
from sgoal import tournament
from sgoal import randbool

############### Generational Genetic Algorithm - GGA ################
# xr: Crossover rate
# xover: Crossover operator
# mutation: Mutation operator
# P: Initial population
# fP: Fitness value for each individual of the population, if provided.
def GGA(f, evals, xr, xover, mutation, P, fP=None):
  n = len(P)
  if( evals>=n and not fP ):
    fP = evaluate_population(f, P)
    evals -= n
  while( evals>=2 ):    
    Q = []
    fQ = []
    for i in range(0,n,2):
      idx1, idx2 = tournament(fP, 2)
      if evals>=2 and randbool(xr):
        a, b = xover(P[idx1], P[idx2])
        a = mutation(a)
        b = mutation(b)
        Q.append(a)
        Q.append(b)
        fQ.append(evaluate(f,a))
        fQ.append(evaluate(f,b))
        evals -= 2
      else:
        Q.append( P[idx1] )
        Q.append( P[idx2] )
        fQ.append(fP[idx1] )
        fQ.append(fP[idx2] )
    P = Q
    fP = fQ
  return P, fP, evals

############### Steady State Genetic Algorithm - GGA ################
# xr: Crossover rate
# xover: Crossover operator
# mutation: Mutation operator
# P: Initial population
# fP: Fitness value for each individual of the population, if provided.
def SSGA(f, evals, xr, xover, mutation, P, fP=None):
  n = len(P)
  if( evals>=n and not fP ):
    fP = evaluate_population(f, P)
    evals -= n
  while( evals>=2 ):    
    if randbool(xr):
      idx1, idx2 = tournament(fP, 2)
      p1, p2 = P[idx1], P[idx2]
      a, b = xover(p1, p2)
      sol = [p1, p2, mutation(a), mutation(b)]
      fsol = [fP[idx1], fP[idx2], evaluate(f,sol[2]), evaluate(f,sol[3])]
      evals -= 2
      fP[idx1], k = best(fsol)
      P[idx1] = sol[k]
      sol.pop(k)
      fsol.pop(k)
      fP[idx2], k = best(fsol)
      P[idx2] = sol[k]      
  return P, fP, evals