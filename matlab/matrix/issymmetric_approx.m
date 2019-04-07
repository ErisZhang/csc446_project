% approximately symmetric up to tolerance `tol`
function b = issymmetric_approx(A, tol)
    if issparse(A)
        A = full(A);
    end
    D = bsxfun(@(x,y) abs(x-y)>tol, tril(A), triu(A)');
    k = find(D);
    b = size(k,1) == 0;
end