% finds spectral radius 
%       \rho(A) = { abs(\lambda) | \lambda is eigenvalue of A }
function r = spectrad(A)
    r = eigs(A,1,'largestabs');
end