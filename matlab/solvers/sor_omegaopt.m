% Gives optimal value for weighted-Jacobi
%   reference: 
%       https://en.wikipedia.org/wiki/Successive_over-relaxation#Convergence
%       Iterative method for Sparse Linear Systems 2ed (p148)
%   
%   Note: not really applicable, has constraints that is not true for stiffness matrix
function opt = sor_omegaopt(A)
    Dinv = diag(diag(A).^-1);
    DinvA = Dinv*A;
    C = sparse(eye(size(A))) - DinvA;
    mu = spectrad(C);
    opt = 1 + (mu / (1 + sqrt(1-mu^2)))^2;
end