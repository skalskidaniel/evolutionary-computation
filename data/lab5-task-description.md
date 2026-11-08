The use of move evaluations (deltas) from previous iterations in local
search
The goal of the task is to improve the time efficiency of the steepest local search with the use move
evaluations (deltas) from previous iterations (list of improving moves) using the neighborhood, which
turned out to be the best in assignment 3. Both inter-route and intra-route moves should be included
in the list. In the case of inter-route moves of the exchange of two edges, you should carefully read
the description of the lectures on the traveling salesman problem, in particular:
• Consider 3 situations (when we browse moves from LM):
• Removed edges (defining the saved move) no longer exist in the current solution (at
least one of them)
• -> remove the move from LM
• Removed edges occur in the current solution in a different relative direction from the
saved one – not applicable now but the move can be applied in the future
• -> leave the move in LM but do not apply it browse LM further
• Removed edges appear in the current solution in the same relative direction (also
both reversed)
• -> perform (apply) the move and remove from LM
• When evaluating new moves we need to consider also moves with relatively inverted edges.
This includes also the initial solution for which we need to evaluate all currently applicable
moves and moves with relatively inverted edges.
This mechanism should be used separately from candidate moves. Optionally, you can try to
implement them both together.
It is advised to use the alternative algorithm using the list of improving moves:
Initiate LM – a list of moves that bring improvement ordered from the best to the worst
Generate an initial solution x
Evaluate all moves applicable to x and add improving moves to LM
repeat
for moves m from LM starting from the best until a applicable move is found
Check if m is applicable and if not remove it from LM
if move m has been found then
x := m(x) (accept m(x))
else
Evaluate all moves applicable to x and add improving moves to LM
until LM is empty
As starting solutions use random solutions.
As baseline report also results of the steepest local search with random starting solutions without
these mechanisms.
Computational experiment: Run each of the methods 200 times.
Reporting results: Use tables as in the previous assignment.
The outline of the report as previously.