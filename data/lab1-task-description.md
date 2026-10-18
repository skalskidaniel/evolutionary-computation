Greedy heuristics
Implement three methods:
• Random solution
• Nearest neighbor considering adding the node only at the end of the current path
• Nearest neighbor considering adding the node at all possible position, i.e. at the end, at the
beginning, or at any place inside the current path
• Greedy cycle
adapted to our problem.
“Nearest” neighbor should not be treated literally. We should always consider the best change of the
objective function, including distance and cost.
For each greedy method generate 200 solutions starting from each node. Generate also 200 random
solutions.
The report should contain:
• Description of the problem
• Pseudocode of all implemented algorithms
• Results of a computational experiment: for each instance and method min, max and average
value of the objective function.
• 2D visualization of the best solution for each instance and method. Cost of nodes should be
presented e.g. by a color, greyscale, or size.
• The best solutions for each instance and method presented as a list of nodes indices (starting
from 0).
• Information whether the best solutions have been checked with the solution checker.
• (Link to) the source code
• Conclusions