


## Problems 


+ global stiffness `K` is neither diagonally dominant nor has a spectral radius of <1 ...
+ how to do gaussian quadrature for matrix integrand?
+ how to correctly enforce dirichlet boundary 
    + set large value to `K_ii` for vertex `i` on boundary
    + set `K_ii` to 1 and `b_i` to 0 for  `i` on boundary, also removes entries `K_{ij},K_{ji}` for all `j`
    + lagrange multiplier
    + projection matrix to contain only unconstrained unknowns
+ convergence speed is slow for classical iterative methods
    + slow for large `A`
        + e^k <= \rho(A) ~ 1 - c/d^2 (nearly 1 for large d) for some constant `c` and `d` is dimension of matrix `A`
+ if cannot get preconditioning to work, work on small


what kind of experiments
+ linear solvers
    + backslash (ground truth)
    + gauss-seidel
    + 


## Meetings 



#### Questions 

+ what are basis function look like in 3d
    + simply 1d hat function in each spatial coordinate ? 
    + complications with using voxels (i.e. box/hexahedral)
        + virtual element methods ?
+ fem enforces natural boundary,
    + imagine 3d mesh, want to find displacement field
    + `u` fixed on feet, but free otherwise 
    + however fdm requires boundary known for the entire boundary, which seem to be impossible to do with fdm?
    + answer: always natural boundary with displacement at nodes not touching ground
        + one more unknown with one more equation 
        + taylor approximation with 3 poitns near boundary to approximate differnetail 
+ what are useful iterative methods (jacobi, conjugate gradient) to use
    + consider parallelism
    + where to look for resources
+ alternative idea
    + if 3d basis for polyhedra is hard to implement, 
        + why not just destruct each voxel into 6 tetrahedra, the local stiffness matrix would be same, 
        + simply need to compute local stiffness matrix 6 times, and populate sparse rhs
        + otherwise, matrix-free iterative method
+ time-wise possible?
    + 3d impl of fem?
    + multigrid?
    + gpu implementation?
    + just pick one!
+ boundary of voxel, model changes ? 
    + need more research on solving pde on voxels 
    + what happens at boundary

#### Meetings

+ fast real=time
    + reconstruct part of the mesh that changes during deformation 
    + initial guess to linear elasticity 
+ linear solver
    + preconditioning: reduce number of iteration for conjugate gradient
        + incomplete cholesky for conjugate gradient solver
    + multigrid methods
        + O(N)
        + might be hard to implement ourselves
+ hex
    + mesh refinement in areas where error is large
        + monitor function, hints to where place are hard
        + i.e. `max_[x_k,x_{k+1}] |x(x)''| `
        + adaptive method for PDE ...
    + adaptive mesh refinement to reduce error ...
        +  i.e. prespecify the grid, or dynamically adjust things itself
+ our idea 
    + solution to previous problem
        + recomputing parts of A
        + supply solution from previous iteration as starting guess ...
    + a better iterative method
        + conjugate gradient 
        + multigrid
+ 446 project
    + error analysis
+ how quickly conjugate gradient converge
    + depends on condition number kappa
    + iteration: sqrt(kappa(A))
    + with preconditioner: sqrt(kappa(M-1 A))
    + monitor residule `Ae = r`, computable with `r=Ax-b`
    + need to find softwares ...




### Junk


```
% per element stress field
se = zeros(12,1);
straine = zeros(size(Tet,1),6);
NV = zeros(size(V,1),1);    % #Vx1 volume of sum of neighboring tets

for i = 1:size(Tet,1)
    B = squeeze(Bs(i,:,:));
    Teti = Tet(i,:);
    ij2p(1:3:end) = 3*Teti-2;
    ij2p(2:3:end) = 3*Teti-1;
    ij2p(3:3:end) = 3*Teti;
    straine(i,:) = B*u(ij2p);
    NV(Teti) = NV(Teti) + Vs(i);
end

% per vertex strain/stress
strain = zeros(size(V,1),6);
stress = zeros(size(V,1),6);

for i =1:size(Tet,1)
    Teti = Tet(i,:);
    strain(Teti,:) = strain(Teti,:) + (1./NV(Teti))*(Vs(i).*straine(i,:));
end

for j = 1:size(V,1)
    stress(j,:) = (C*strain(j,:)')';
end

% von mises
se = zeros(6,1);
VM = zeros(size(V,1),1);
for j = 1:size(V,1)
    se = stress(j,:);
    VM(j) = (1./sqrt(2))*sqrt((se(1)-se(2))^2 + (se(2)-se(3))^2 + ...
        (se(3)-se(1))^2 + 6*(se(4)^2+se(5)^2+se(6)^2));
end
```