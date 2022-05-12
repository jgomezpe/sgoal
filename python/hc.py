from sgoal import pick
from sgoal import evaluate

# Classical Hill Climbing Algorithm with neutral mutations
# f: Function to be optimized
# variation: Variation operation
# evals: Maximum number of fitness evaluations
# x: Initial point
# fx: f value at point x (if provided)
def HC( f, variation, evals, x, fx=None):
  if(evals>0 and not fx): 
    fx = evaluate(f, x)
    evals-=1
  for i in range(evals):
    y = variation(x)
    fy = evaluate(f,y)
    x, y, fx, fy = pick( x, y, fx, fy )
  return x, fx
