# Newton in RefNewton
Newton in RefNewton does not converge for large Jacobians.

# Random Walk
Random walk is not walking (working).

# Subgrid-Scale
Sometimes gives infinite pushes. Also lacks proper high-order accurate time integration (Anna, Marcel say it might be okay).

# RefMapping
Random crashes on Hawk when running huge grids.<br/>
```
Tolerance Issue with internal element
halo-elem = T
```

# ParaView
## Parallel VTU
Distributed Data Filter does not correctly find the global node IDs and crashes due to missing points. We need to build the unique global node IDs and give them to ParaView. Current idea is to modify CreateConnectivity, so MPI root builds the global IDs and distributes them to the procs.

## Surface Visu
Currently, SurfStates can only be visualized with paraviewposti.<br/>
- Step 1: Move SurfState to normal state file.<br/>
- Step 2: Create a multiblock reader for the supersampled surface faces.<br/>
