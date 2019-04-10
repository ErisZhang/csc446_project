% weighted Jacobi Method for {s}parse matrix
%
%   Convergence;
%       A SPD && 0 < \oemga < 2/\lamba_{max}(D^{-1}A)
%       https://en.wikipedia.org/wiki/Jacobi_method
%
%   Iterations:
%       x^{k+1}_{i} = \omega * (1/a_{ii}) [ \sum_{j=1,j\neq i}^n -a_{ij}*x_{j}^{k} + b_i ]
%                   + (1-\Omega) * x^{k}_{i}
%
function x = jacobis(A,b,x0,max_iter,tol,omega)
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
    Dinv = D.^-1;

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
            x(i) = omega * Dinv(i) * ( -A(i,nnzi)*xp(nnzi) + A(i,i)*xp(i) + b(i) ) + ...
               (1-omega) * xp(i);
        end

        r = (x-xp).*Dinv;

        if norm(r,inf) <= tol
            return;
        end
    end
end