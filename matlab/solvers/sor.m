% Successive-over relaxation on {s}parse matrix
%
%       x_{i}^{k+1} = (omega/a_{ii}) * [
%           -\sum_{j=1}^{i-1} a_{ij} x_{j}^{k+1} - \sum_{j=i+1}^n a_{ij} x_{j}^{k} + b_{i}]
%           + (1-omega) * x_{i}^{k}
%
function x = sor(A,b,x0,max_iter,tol,omega)
    if  ~issparse(A)
        A = sparse(A);
        warning('A should be sparse');
    end

    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');
    assert(all(full(diag(A)~=0)), 'A_{ii} cannot be 0 for all i');

    nnzibyrow = nnz_indices_byrow(A);

    n = size(A,1);
    D = diag(A);
    Dinv = full(D.^-1);

    nnzi = zeros(size(nnzibyrow{1}));
    r  = zeros(size(x0));
    x  = x0;
    xp = x0;

    for it = 1:max_iter
        xp = x;

        if mod(it, 10) == 0
            [it max_iter norm(r,inf)]
        end

        for i = 1:n
            nnzi = nnzibyrow{i};
            x(i) = omega * Dinv(i) * ( -A(i,nnzi)*x(nnzi) + A(i,i)*x(i) + b(i) ) + ...
               (1-omega) * xp(i);
        end

        % Compute residual
        %   https://www3.nd.edu/~zxu2/acms40390F11/sec7-3-2.pdf
        r = (x-xp).*Dinv;

        if norm(r,inf) <= tol
            return;
        end
    end
end