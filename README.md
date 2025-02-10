# sgoal
A global optimization problem is formulated in terms of finding a point x in a subset Ω ⊆ Φ where a certain function f : Φ → R, attains is best/optimal value (minimum or maximum). In the optimization field, Ω, Φ, and f are called the feasible region, the solution space, and the objective function, respectively. The optimal value for the objective function (denoted as f<sup>*</sup> ∈ R) is
suppose to exist and it is unique (R is a total order).

A Stochastic Global Optimization ALgorithm (SGoal) is an iterative algorithm that generates a new candidate set of solutions (called population) from a given population using a stochastic operation 

A formal description is provided by Jonatan Gomez in "Stochastic global optimization algorithms: A systematic formal approach", Information Sciences, Volume 472, 2019, Pages 53-76, ISSN 0020-0255, https://doi.org/10.1016/j.ins.2018.09.021. (https://www.sciencedirect.com/science/article/pii/S0020025517305248)

This library includes the following SGoals:
<ul>
  <li>HC: Hill Climbing with neutral mutations.</li>
  <li>GGA: Generational Genetic Algorithm with tournament 4 as parent's selection mechanism.</li>
  <li>SSGA: Steady State Genetic Algorithm with tournament 4 as parent's selection mechanism.</li>
  <li>CHAVELA: Canonical Hibrid Adaptive EVolutionary ALgorithm.</li>
</ul>
