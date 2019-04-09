% Finds boundary vertices s.t.
%       V(i,dim) < ratio*max_length_along(dim)
function b = boundary_vertices(V,dim,ratio)

    b = zeros(0,1);
    nbc = 1;
    
    maxd = max(V(:,dim));
    mind = min(V(:,dim));
    tol = ratio * abs(maxd-mind);

    for i = 1:size(V,1)
        if abs(V(i,dim)-mind) < tol
            b(nbc) = i;
            nbc = nbc + 1;
        end
    end
end