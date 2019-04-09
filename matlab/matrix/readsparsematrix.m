% reads in sparse matrix in a txt file, of form
%       i   j   v
%       ...
function S = readsparsematrix(filename)
    IJV = readmatrix(filename,'Delimiter',' ');
    assert(size(IJV,2) == 3, 'readsparsematrix: ncols must be 3');
    S = sparse(IJV(:,1),IJV(:,2),IJV(:,3));
end