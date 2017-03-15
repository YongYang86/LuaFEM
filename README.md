# LuaFEM

##Highlights
poisson_p1.lua

 1. Read P1 or Q1 finite element from .msh file
 2. Save mesh as a .vtk file
 3. Assemble the mass matrix
 4. Assemble the lumped mass matrix
 5. Assemble the stiffness matrix $ (\nabla\phi_i,\nabla\phi_j)$
 6. CG iterative method for sparse matrix Ax=b
##Improvements
 1. Create structured mesh without reading .msh file.
 2. GMRES iterative method

##Problems
 1. If one uses "luajit possion_p1.lua", it becomes fast. However, one meets the problem "luajit: not enough memory".
## Author

Yong Yang (wacyyang@gmail.com).

## License

See LICENSE file.

