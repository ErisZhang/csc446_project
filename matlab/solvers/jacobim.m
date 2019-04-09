% Jacobi iterative method for `Ax = b` using {m}atrix multiply
%       splitting:
%           A = D - R       where `R = -(tril(A)+triu(A))`
%       iteration:      x^{k+1} = D^{-1} (R*x^{k} + b)
%                               = H*x^{k} + v
%           H = D^{-1}*R
%           v = D^{-1}*b
%
%       convergence: \rho(A) < 1 && diagonally dominant
%       reference: https://ece.uwaterloo.ca/~dwharder/NumericalAnalysis/04LinearAlgebra/jacobi/
%
function x = jacobim(A,b,x0,max_iter,tol)
    assert(isdiagdominbyrow(A), 'A needs to be diagonally dominant');
    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');

    D = diag(A);
    R = diag(D) - A;

    Dinv = diag(D.^-1);
    H = Dinv*R;
    v = Dinv*b;

    x  = x0;
    xp = x0;

    for i = 1:max_iter
        xp = x;
        x = H*x + v;
        if norm(x-xp,inf) < tol
            return;
        end
    end
    warning('jacobi iteration did not converge');
end