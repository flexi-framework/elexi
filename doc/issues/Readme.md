# Newton in RefNewton
Newton in RefNewton does not converge for large Jacobians.

# Random Walk
Random walk is not walking (working).

# Subgrid-Scale
Sometimes gives infinite pushes. Also lacks proper high-order accurate time integration (Anna, Marcel say it might be okay).

# ParaView
## Parallel VTU
Parallel VTUs can only be read with a number of processors smaller than the number of corresponding .vtu files. Otherwise, random gradients appear messed up.
## Surface Visu
Currently, SurfStates can only be visualized with paraviewposti.<br/>
- Step 1: Move SurfState to normal state file.<br/>
- Step 2: Create a multiblock reader for the supersampled surface faces.<br/>
