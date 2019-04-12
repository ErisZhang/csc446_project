% Successive-over relaxation on {s}parse matrix
%
%       Convergence: A SPD 
%           reference: Iterative method for Sparse Linear Systems 2ed (p105)
%
%       x_{i}^{k+1} = (omega/a_{ii}) * [
%           -\sum_{j=1}^{i-1} a_{ij} x_{j}^{k+1} - \sum_{j=i+1}^n a_{ij} x_{j}^{k} + b_{i}]
%           + (1-omega) * x_{i}^{k}
%
%   Inputs
%       varargin:
%           'Omega', (defaults to 1)
%                    0 < omega < 2
%
%   Outputs
%       varargout:
%           xks     each coloum is x^k
%           rks     each coloum is r^k
%
function [x,varargout] = sor(A,b,x0,max_iter,tol,varargin)
    assert(issparse(A), 'A has to be sparse');
    assert(numel(unique([size(A,2) size(b,1) size(x0,1)])) == 1, ...
        'A,b dimension should match');

    % Gauss-Seidel iteration
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

    [x,xks,rks] = sors(A,b,x0,max_iter,tol,omega,saveon);
    varargout{1} = xks;
    varargout{2} = rks;
end




% SOR on sparse `A`
function [x,varargout] = sors(A,b,x0,max_iter,tol,omega,saveon)

    nnzibyrow = nnz_indices_byrow(A);

    n = size(A,1);
    D = diag(A);
    Dinv = full(D.^-1);

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
            fprintf('(SOR) %d/%d:\t%.5e (rel residual)\n', it,max_iter,norm(r)/normb);
        end

        for i = 1:n
            nnzi = nnzibyrow{i};
            x(i) = omega * Dinv(i) * ( -A(i,nnzi)*x(nnzi) + A(i,i)*x(i) + b(i) ) + ...
                (1-omega) * xp(i);
        end

        r = b - A*x;

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


% Following code not really in use ..


% Gauss-Seidel iterative method for `Ax = b` where `A` is {s}parse
%
%       x_{i}^{k+1} = (1/a_{ii}) * [
%           -\sum_{j=1}^{i-1} a_{ij} x_{j}^{k+1} - \sum_{j=i+1}^n a_{ij} x_{j}^{k} + b_{i}
%       ];
%
%       Compute residual
%           https://www3.nd.edu/~zxu2/acms40390F11/sec7-3-2.pdf
%       Relative residual 
%           https://www8.cs.umu.se/kurser/5DA002/HT08/reassignment4.pdf
%
function x = gaussseidel(A,b,x0,max_iter,tol)
    if issparse(A)
        x = gaussseidels(A,b,x0,max_iter,tol);
    else
        x = gaussseideld(A,b,x0,max_iter,tol);
    end
end


% Gauss-Seidel iterative method on sparse `A`
function x = gaussseidels(A,b,x0,max_iter,tol)

    nnzibyrow = nnz_indices_byrow(A);

    n = size(A,1);
    D = diag(A);
    Dinv = full(D.^-1);
    
    nnzi = zeros(size(nnzibyrow{1}));
    r  = zeros(size(x0));
    x  = x0;
    xp = x0;
    normb = norm(b);

    for it = 1:max_iter
        xp = x;

        for i = 1:n
            nnzi = nnzibyrow{i};
            x(i) = Dinv(i) * ( -A(i,nnzi)*x(nnzi) + A(i,i)*x(i) + b(i) );
        end

        r = b - A*x;

        if norm(r)/normb <= tol
            return;
        end
    end
end


% Gauss-Seidel iterative method on dense `A`
function x = gaussseideld(A,b,x0,max_iter,tol)

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

        r = b - A*x;
    
        if norm(r)/normb <= tol
            return;
        end
    end
end