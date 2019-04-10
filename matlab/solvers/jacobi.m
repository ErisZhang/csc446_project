% weighted Jacobi Method
%
%   Convergence;
%       A SPD && 0 < \oemga < 2/\lamba_{max}(D^{-1}A)
%       https://en.wikipedia.org/wiki/Jacobi_method
%
%   Iterations:
%       x^{k+1}_{i} = \omega * (1/a_{ii}) [ \sum_{j=1,j\neq i}^n -a_{ij}*x_{j}^{k} + b_i ]
%                   + (1-\Omega) * x^{k}_{i}
%   Inputs
%       varargin:
%           'Omega', (defaults to 1)
%                    0 < omega < 1
%
function x = jacobi(A,b,x0,max_iter,tol,varargin)
    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');

    % un-weighted jacobi
    omega = 1;

    % Map of parameter names to variable names
    params_to_variables = containers.Map( ...
        {'Omega'}, ...
        {'omega'});
    v = 1;
    while v <= numel(varargin)
        param_name = varargin{v};
        if isKey(params_to_variables,param_name)
            assert(v+1<=numel(varargin));
            v = v+1;
            % Trick: use feval on anonymous function to use assignin to this workspace
            feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
        else
            error('Unsupported parameter: %s',varargin{v});
        end
        v=v+1;
    end

    if issparse(A)
        x = jacobis(A,b,x0,max_iter,tol,omega);
    else
        x = jacobid(A,b,x0,max_iter,tol);
    end
end


% Weighted-Jacobi on sparse `A`
function x = jacobis(A,b,x0,max_iter,tol,omega)

    nnzibyrow = nnz_indices_byrow(A);

    n = size(A,1);
    D = diag(A);
    Dinv = D.^-1;

    nnzi = zeros(size(nnzibyrow{1}));
    r  = zeros(size(x0));
    x  = x0;
    xp = x0;
    normb = norm(b);

    for it = 1:max_iter
        xp = x;

        % if mod(it, 10) == 0
        %     [it max_iter norm(r,inf)]
        % end

        for i = 1:n
            nnzi = nnzibyrow{i};
            x(i) = omega * Dinv(i) * ( -A(i,nnzi)*xp(nnzi) + A(i,i)*xp(i) + b(i) ) + ...
               (1-omega) * xp(i);
        end

        r = (x-xp).*Dinv;

        if norm(r)/normb <= tol
            return;
        end
    end
end


% Jacobi iterative method on dense `A`
%       iteration:
%           x^{k+1}_{i} = (1/a_{ii}) [ \sum_{j=1,j\neq i}^n -a_{ij}*x_{j}^{k} + b_i ]
%
function x = jacobid(A,b,x0,max_iter,tol)
    assert(all(full(diag(A)~=0)), 'A_{ii} cannot be 0 for all i for 1-Jacobi');

    n = size(A,1);
    D = diag(A);
    Dinv = D.^-1;

    r  = zeros(size(x0));
    x  = x0;
    xp = x0;
    normb = norm(b);

    for it = 1:max_iter
        xp = x;

        for i = 1:n
            x(i) = Dinv(i) * ( -A(i,:)*x + A(i,i)*x(i) + b(i) );
        end

        r = (x-xp).*Dinv;
    
        if norm(r)/normb <= tol
            return;
        end
    end
end