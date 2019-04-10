% Returns cell `C` for sparse matrix `S` where
%       S(i,C{i}(j)) is `j`-th nonzero value in row `i`
%   Reference:
%   https://www.mathworks.com/matlabcentral/answers/344200-find-out-non-zero-indexes-of-rows-in-matrix
%
function indices = nnz_indices_byrow(S)
    assert(issparse(S), 'S has to be sparse');
    [c, r] = find(S.');
    indices = accumarray(r, c, [size(S,1), 1], @(L) {L.'} );
end 