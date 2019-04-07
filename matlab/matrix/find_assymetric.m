% Finds (i,j) s.t. A_{ij} ~= A_{ji} && i > j (lower triangle)
%       idx     #nnz x 2
function IJ = find_assymetric(A)
    if issymmetric(A)
        idx = []; return;
    end

    if issparse(A)
        A = full(A);
    end

    D = bsxfun(@(x,y) x~=y, tril(A), triu(A)');
    k = find(D);
    [I,J] = ind2sub(size(A), k);

    IJ = [I,J];
end