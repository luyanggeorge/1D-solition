# Single-Line Solition Test

- 2D simulation on single-line soliton through PFE. The system is initialised with exact solution of KPE for single solition.
- This version has not implement periodic boundary condition yet.
- Before running the simulaton, `KPE_solution.py` can be used to show the profile of the wave at the start and the end time steps.
- Codes start with `pp` are for post-processing of numerical results.

## Comparison between two ways to initialise the velocity potentials
Strictly speaking, we should use the z coordinates before coordinate transformation when initialising the velocity potentials. For example, for the velocity potential at the free surface, z should be substituted with h rather than H0. However, It can be seen from the results that this doesn't afftect the performance. 
