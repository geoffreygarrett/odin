Typically, in the context of Monte Carlo Tree Search (MCTS), every visit to a node involves an
update of that node's statistics. Therefore, most MCTS implementations keep a single count that
serves both as the number of visits to a node and the number of updates to that node's statistics.

Here's a high-level explanation of why this is the case:

- **Selection:** Starting from the root, the tree is traversed to a leaf node. Each node along the
  path increments its visit count.

- **Expansion:** The leaf node is expanded (if possible) by adding one or more child nodes. The
  child node (or nodes) represents a state not previously explored from the current state.

- **Simulation:** From the newly added child node, a simulation (or "rollout") is run to a terminal
  state. The result of this simulation represents an estimate of the value of the new state.

- **Backpropagation:** The result of the simulation is "backpropagated" up the tree, updating the
  statistics of each node along the path from the new child node to the root.

During the backpropagation phase, each node in the path from the newly added node up to the root has
its statistics updated. This includes incrementing its visit count and updating its value estimate (
usually, the average of the simulation results for each visit).

So, in conclusion, in standard MCTS, there's no need to keep separate visit and update counts
because every visit results in an update. However, depending on your specific use case or any
modifications to the standard MCTS algorithm, you might need to keep track of these separately.