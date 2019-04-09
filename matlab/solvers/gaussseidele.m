% Gauss-Seidel iterative method for `Ax = b` on each {e}lement
%
%       x_{i}^{k+1} = (1/a_{ii}) * [
%           -\sum_{j=1}^{i-1} a_{ij} x_{j}^{k+1} - \sum_{j=i+1}^n a_{ij} x_{j}^{k} + b_{i}
%       ];
%
function x = gaussseidele(A,b,x0,max_iter,tol)
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
            x(i) = Dinv(i) * ( -A(i,:)*x + A(i,i)*x(i) + b(i) );
        end

        if norm(x-xp,inf) < tol
            return;
        end
    end
end


