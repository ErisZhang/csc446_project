% Jacobi iterative method for `Ax = b` on each {e}lement
%       iteration:
%           x^{k+1}_{i} = (1/a_{ii}) [ \sum_{j=1,j\neq i}^n -a_{ij}*x_{j}^{k} + b_i ]
%
function x = jacobie(A,b,x0,max_iter,tol)
    assert(isdiagdominbyrow(A), 'A needs to be diagonally dominant');
    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');

    n = size(A,1);
    D = diag(A);
    Dinv = D.^-1;

    x  = x0;
    xp = x0;

    for it = 1:max_iter
        xp = x;

        for i = 1:n
            x(i) = Dinv(i) * ( -A(i,:)*xp + A(i,i)*xp(i) + b(i) );
        end

        if norm(x-xp,inf) < tol
            return;
        end
    end
end



