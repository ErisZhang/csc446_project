

## Junk 

% Gauss-Seidel iterative method for `Ax = b` using {m}atrix multiply
%       splitting:
%           A = (D-L) - (U)
%       iteration:      x^{k+1} = (D-L)^{-1}*U*x^{k} + (D-L)^{-1}*b;
%
%       convergence: A symmetric positive definite
%
% Note matrix does not work!

% function x = gaussseidelm(A,b,x0,max_iter,tol)
%     assert(isdiagdominbyrow(A), 'A needs to be diagonally dominant');
%     assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
%         'A,b dimension should match');

%     D = diag(diag(A));
%     L = -tril(A);
%     U = -triu(A);

%     Pinv = inv(D-L);
%     H = Pinv*U;
%     v = Pinv*b;

%     x  = x0;
%     xp = x0;

%     for i = 1:max_iter
%         xp = x;
%         x = H*x + v;
%         if norm(x-xp,inf) < tol
%             return;
%         end
%     end
%     warning('gauss-seidel iteration did not converge');
% end




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

    Dinv = D.^-1;
    H = diag(Dinv)*R;
    v = diag(Dinv)*b;

    r  = zeros(size(x0));
    x  = x0;
    xp = x0;

    for i = 1:max_iter
        xp = x;
        x = H*x + v;
        r = (x-xp).*Dinv;
        if norm(r,inf) <= tol
            return;
        end
    end
end

