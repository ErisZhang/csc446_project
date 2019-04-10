% if matrix is diagonally dominant by row
function [b,l] = isdiagdominbyrow(A)
    if issparse(A)
        assert(true)
    else
        absA = abs(A);
        l = 2*diag(absA) > sum(absA,2);
        b = all(l);
    end
end