# Assignment 3: Local Search

## Local search
Implement local search. Local search should use deltas (changes) of the objective function. Do not use trivial implementation with neighbor solution constructed explicitly and the objective function calculated from scratch. On the other hand, do not use more advanced techniques like deltas from previous iterations or candidate moves (we will implement them in future assignments).

## Type of local search
Implement both steepest and greedy version of local search. In greedy version the neighborhood should be browsed in random/randomized order. In the report please describe what kind of randomization was used. Note, that in greedy version, it is not acceptable to evaluate (calculate delta) all moves and them browse them in random order and select the first improving one, because we need to minimize the number of evaluated moves.

## Type of neighborhood
We will need two kinds of moves: 
* **intra-route moves** - moves changing the order of nodes within the same set of selected nodes.
* **inter-route moves** - moves changing the set of selected nodes. 

As inter-route move use exchange of two nodes one selected one not selected. 

For intra-route moves use two options:
* two-nodes exchange,
* two-edges exchange.

As inter-route moves use always nodes-exchange with one selected one not selected.

Note that the whole neighborhood is composed of moves of two kinds inter- and intra-route moves. In the steepest version we should select the best move from the whole neighborhood composed of moves of two types. In the greedy version ideally we should test moves of two kinds in a random order. If you use a simpler randomization do not, for example, start always from moves of a given kind (e.g. always from intra-route moves).

## Type of starting solutions
Use two options:
* Random starting solutions.
* The best greedy construction heuristic (including regret heuristics) from the previous assignments with a random starting node.

Summarizing we have 8 combinations based on three binary options: type of local search (steepest or greedy), type of intra-route moves (nodes or edges exchange), type of starting solutions (random or greedy).

## Computational experiment
Run each of the eight methods 200 times. For random starting solutions use 200 randomly generated solutions. For greedy starting solutions use each of the 200 nodes as the starting node for the greedy heuristic.

## Reporting results
Use table with the following outline:

| | Instance 1 | Instance 2 |
| :--- | :--- | :--- |
| **Method 1** | av (min-max) | av (min-max) |
| **Method 2** | av (min - max) | av (min-max) |

Include also results for all previous methods.

Report two tables one for the values of the objective function and for the running times. 

The outline of the report as previously.


extra guidelines

Do not use methods more advanced than delta

the evaluations shall show that Greedy is faster than steepest. If not the method is implemented badly

Describe how we perform the randomization in the report

Pseudocode for local search should be precise, others very brief (summarized)

Also add the calculation time of the initial solution

Improve the visualizations- no red dots at the visited points, the line connecting them shall be black, the points should have the colour indicating their price