# sgoal
A global optimization problem is formulated in terms of finding a point x in a subset Ω ⊆ Φ where a certain function f : Φ → R, attains is best/optimal value (minimum or maximum). In the optimization field, Ω, Φ, and f are called the feasible region, the solution space, and the objective function, respectively. The optimal value for the objective function (denoted as f<sup>*</sup> ∈ R) is
suppose to exist and it is unique (R is a total order).

A Stochastic Global Optimization ALgorithm (SGoal) is an iterative algorithm that generates a new candidate set of solutions (called population) from a given population using a stochastic operation 

A formal description is provided by Jonatan Gomez in "Stochastic global optimization algorithms: A systematic formal approach", Information Sciences, Volume 472, 2019, Pages 53-76, ISSN 0020-0255, https://doi.org/10.1016/j.ins.2018.09.021. (https://www.sciencedirect.com/science/article/pii/S0020025517305248)

This library includes the following SGoals:
<ul>
  <li>GGA: Generational Genetic Algorithm with tournament 4 as parent's selection mechanism.</li>
  <li>SSGA: Steady State Genetic Algorithm with tournament 4 as parent's selection mechanism.</li>
  <li>CHAVELA: The Canonical Hibrid Adaptive EVolutionary ALgorithm as proposed by Jonatan Gómez and Elizabeth León in: "On the class of hybrid adaptive evolutionary algorithms (chavela)" Published in Natural Computing / Issue 3/2021, Print ISSN: 1567-7818 Electronic ISSN: 1572-9796, DOI (https://doi.org/10.1007/s11047-021-09843-5)</li>
  <li>HC: Hill Climbing with neutral mutations.</li>
  <li>Rule_1_5th: 1+1 Evolutionary Strategy (Hill Climbing) with neutral mutations and 1/5th rule, see Beyer, Hans-Georg & Schwefel, Hans-Paul. (2002). Evolution strategies - A comprehensive introduction. Natural Computing. 1. 3-52. 10.1023/A:1015059928466.</li>
  <li>RMHC: The HC algorithm suggested by Richard Palmer, that Forrest and Mitchell named as "random mutation hill-climbing" (RMHC), see M. Mitchell and J. Holland, “When will a genetic algorithm outperform hill-climbing?”. Santa Fe Institute, Working Papers, 01 1993.</li>
  <li>GGS1: The Global Search Algorithm for GA-easy functions propossed by Das and Whitley 1991, which tries only order 1-schemas, see R. Das and L. D. Whitley, “The only challenging problems are deceptive: Global search by solving order-1 hyperplanes, in Proceedings of the 4th International Conference on Genetic Algorithms, San Diego, CA, USA, July 1991 (R. K. Belew and L. B. Booker, eds.), pp. 166– 173, Morgan Kaufmann, 1991</li>
  <li>GGSC1: A Generalization of The Global Search Algorithm with complement propossed by G. Venturini 1995 in First International Conference on Genetic Algorithms in Engineering Systems: Innovations and Applications, pp. 357–364, 1995.</li>
  <li>GABO: The "GABO: Gene Analysis Base Optimization" algorithm as proposed by J. Gomez and E. Leon, "Gabo: Gene Analysis Bitstring Optimization," In 2022 IEEE Congress on Evolutionary Computation (CEC), Padua, Italy, 2022, pp. 1-8, doi: 10.1109/CEC55065.2022.9870237.</li>
</ul>
