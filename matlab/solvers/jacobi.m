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
%           'SaveOn', []
%                   indices `k`s' for which {x^k, r^k} will be saved
%
%   Outputs
%       varargout:
%           xks     each coloum is x^k
%           rks     each coloum is r^k
%
function [x,varargout] = jacobi(A,b,x0,max_iter,tol,varargin)
    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');

    % un-weighted jacobi
    omega = 1;
    saveon = zeros(0,1);

    % Map of parameter names to variable names
    params_to_variables = containers.Map( ...
        {'Omega','SaveOn'}, ...
        {'omega','saveon'});
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
        [x,xks,rks] = jacobis(A,b,x0,max_iter,tol,omega,saveon);
    else
        [x,xks,rks] = jacobid(A,b,x0,max_iter,tol,saveon);
    end
    varargout{1} = xks;
    varargout{2} = rks;
end


% Weighted-Jacobi on sparse `A`
%       omega = 2/3 is a good value
function [x,varargout] = jacobis(A,b,x0,max_iter,tol,omega,saveon)

    nnzibyrow = nnz_indices_byrow(A);

    n = size(A,1);
    D = diag(A);
    Dinv = D.^-1;

    nnzi = zeros(size(nnzibyrow{1}));
    r  = zeros(size(x0));
    x  = x0;
    xp = x0;
    normb = norm(b);

    saveoni = 1;
    xks = zeros(size(x0,1),0);
    rks = zeros(size(x0,1),0);

    for it = 1:max_iter
        xp = x;

        if mod(it, 50) == 0
            fprintf('(Jacobi) %d/%d:\t%.5e (rel residual)\n',it,max_iter,norm(r)/normb);
        end

        for i = 1:n
            nnzi = nnzibyrow{i};
            x(i) = omega * Dinv(i) * ( -A(i,nnzi)*xp(nnzi) + A(i,i)*xp(i) + b(i) ) + ...
               (1-omega) * xp(i);
        end

        r = (x-xp).*Dinv;

        if ~isempty(saveon) && saveoni <= max(size(saveon)) && saveon(saveoni) == it
            xks(:,saveoni) = x;
            rks(:,saveoni) = r;
            saveoni = saveoni + 1;
        end

        if norm(r)/normb <= tol
            break;
        end
    end

    varargout{1} = xks;
    varargout{2} = rks;
end


% Jacobi iterative method on dense `A`
%       iteration:
%           x^{k+1}_{i} = (1/a_{ii}) [ \sum_{j=1,j\neq i}^n -a_{ij}*x_{j}^{k} + b_i ]
%
function [x,varargout] = jacobid(A,b,x0,max_iter,tol,saveon)
    assert(all(full(diag(A)~=0)), 'A_{ii} cannot be 0 for all i for 1-Jacobi');

    n = size(A,1);
    D = diag(A);
    Dinv = D.^-1;

    r  = zeros(size(x0));
    x  = x0;
    xp = x0;
    normb = norm(b);

    saveoni = 1;
    xks = zeros(size(x0,1),0);
    rks = zeros(size(x0,1),0);

    for it = 1:max_iter
        xp = x;

        for i = 1:n
            x(i) = Dinv(i) * ( -A(i,:)*x + A(i,i)*x(i) + b(i) );
        end

        r = (x-xp).*Dinv;

        if ~isempty(saveon) && saveoni <= max(size(saveon)) && saveon(saveoni) == it
            xks(:,saveoni) = x;
            rks(:,saveoni) = r;
            saveoni = saveoni + 1;
        end
    
        if norm(r)/normb <= tol
            break;
        end
    end

    varargout{1} = xks;
    varargout{2} = rks;
end