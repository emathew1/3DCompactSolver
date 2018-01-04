# 3DCompactSolver
A compressible Navier-Stokes solver discretized using a 6th order compact method in space (Lele, 1992) and RK4 in time written in C++. This version is only for use on shared-memory machines as of now, multi-threading implemented via OpenMP.

Current implementation assumes ideal gas, diffusive terms have been expanded to take advantage of the resolution properties of the discretization method, and an 8th order compact filter is to remove the high wavenumber instabilities due to the central discretization while retaining the order of accuracy of the method.  

Currently only periodic, no-slip adiabatic wall, moving adiabatic wall, and non-reflective sponge boundary conditions have been implemented. Geometries in the solver is restricted to cartesian meshes and uniform X, Y, or Z spacing for now. 
