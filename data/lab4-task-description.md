# The Use of Candidate Moves in Local Search

## Goal
The objective of this task is to improve the time efficiency of the **steepest local search** by using candidate moves. You will use the neighborhood definition that proved to be the best in the previous assignment. 

While the quality of the solutions might be slightly worse when using candidate moves, the execution time should significantly decrease.

## Definitions
* **Candidate Moves:** Moves that introduce at least one candidate edge to the solution. 
* **Candidate Edges:** For each vertex, define its 10 "nearest" vertices. "Nearest" is calculated by taking into account the sum of the edge length and the vertex cost. 
  * *Note:* The parameter for the number of nearest vertices (defaulted to 10) can be selected experimentally to obtain the best results.

## Move Types
Both move types must ensure that at least one candidate edge is introduced to the solution. Note that in both cases, there are two possible moves that introduce a given candidate edge.

* **Intra-route move:** Use **2-opt** (e.g., two-edges exchange).
* **Inter-route move:** Introduce a new move consisting of inserting a new (unselected) vertex next to a selected, existing one (node exchange). 
  * *Crucial constraint:* It is not correct to exchange a selected node with one of its nearest neighbors without adding at least one candidate edge to the solution.

## Implementation Details
Since the goal is to improve time efficiency, **it is not acceptable to iterate over all possible moves** and check if a move is a candidate or not. You must iterate only over candidate moves. 

**Suggested Loop Structure:**
1. Create a first loop iterating over all selected nodes.
2. Create a second inner loop iterating over the nearest candidate nodes to the node from the first loop.
3. Depending on whether the candidate node is currently selected or not, evaluate the appropriate candidate moves.

## Experimental Setup
* **Starting Solutions:** Use random solutions.
* **Baseline:** Report the results of the standard steepest local search (using random starting solutions) *without* the candidate moves mechanism for comparison.
* **Computational Experiment:** Run each method **200 times**.

## Reporting Results
* Use tables formatted in the same way as the previous assignment.
* Follow the same general outline for the report as previously used.