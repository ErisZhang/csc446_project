% Finds (i,j) s.t. A_{ij} \not\approx A_{ji} && i > j (lower triangle)
%       idx     #nnz x 2
function IJ = find_assymetric_approx(A, tol)
    if issymmetric(A)
        idx = []; return;
    end

    if issparse(A)
        A = full(A);
    end

    D = bsxfun(@(x,y) abs(x-y)>tol, tril(A), triu(A)');
    k = find(D);
    [I,J] = ind2sub(size(A), k);

    IJ = [I,J];
end