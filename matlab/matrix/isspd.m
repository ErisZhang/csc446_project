% checks if nxn matrix `A` is symmetric positive definite
function b = isspd(A)
    b = ispd(A) && issymmetric(A);
end