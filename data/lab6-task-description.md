Multiple start local search (MSLS) and iterated local search (ILS)
The goal of the task is to implement two simple extensions of local search:
• Multiple start local search (MSLS) – we will use steepest local search starting from random
solutions.
• Iterated local search (ILS).
You shall use the version with list of
moves of steepest local search.
The perturbation for ILS should be designed by you and precisely described in the report. You should
aim, however, at obtaining better results with ILS than with MSLS.
Computational experiment: Run each of the methods (MSLS and ILS) 20 times for each instance. In
MSLS perform 400 iterations of basic local search. For ILS as the stopping condition use the average
running time of MSLS. For ILS report also the number of runs of basic local search. For ILS as the
starting solution (one for each run of ILS) use random solution. Note that the results of a single run of
MSLS is the best solution among a given number of runs of local search.
Reporting results: Use tables as in the previous assignment. For ILS add a table with the number of
runs of basic LS. Report the best solutions for each instance as a lists of nodes. Include results of the
best construction heuristics(s) and the best version(s) of local search.
The outline of the report as previously.
The pesudocode of local search doesnt need to be very detailed, but you shall state clearly what was the perturbation used in the ILS and describe its behaviour.
