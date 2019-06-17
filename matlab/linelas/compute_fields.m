% for computing per-vertex fields `VM` given 3nx1 array `u`
%       either solution or residual to linear system
function [U,strain,stress,vm,VM] = compute_fields(V,Tet,C,Bs,u)
    U = zeros(size(V));
    U(:,1) = u(1:3:end);
    U(:,2) = u(2:3:end);
    U(:,3) = u(3:3:end);

    [strain,stress,vm]=per_element_fields(Tet,C,Bs,u);

    N = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        N(Tet(i,:)) = N(Tet(i,:)) + 1;
    end
    VM = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        VM(Tet(i,:),1) = VM(Tet(i,:),1) + vm(i)./N(Tet(i,:),1);
    end
end

