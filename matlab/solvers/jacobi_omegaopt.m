% Gives optimal value for weighted-Jacobi
%   reference: https://en.wikipedia.org/wiki/Jacobi_method#Weighted_Jacobi_method
%
function opt = jacobi_omegaopt(A)
    Dinv = diag(diag(A).^-1);
    DinvA = Dinv*A;
    eigsmax = eigs(DinvA,1,'largestabs');
    eigsmin = eigs(DinvA,1,'smallestabs');
    opt = 2/(eigsmax+eigsmin);
end