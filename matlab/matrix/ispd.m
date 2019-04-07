% checks if nxn matrix `A` is positive definite
function y = ispd(A)
    [~,flag] = chol(A);
    y = flag == 0;
end