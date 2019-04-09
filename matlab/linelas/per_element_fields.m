% Computes per-element stress/strain/vm fields
%
% Inputs:
%   Tet     #Tet x 4
%   C       12 x 12         stiffness tensor (`\sigma = C * \epsilon`)
%                                   
%   Bs      #Tet x 6 x 12   B matrix for all tets
%   u       #V*3 x 1        node-wise ordered vertex displacement
function [strain,stress,vm] = per_element_fields(Tet, C, Bs, u)
    Teti = zeros(1,4);
    ij2p = zeros(12,1); % for indexing from `V` to `K/f/u`
    B = zeros(6,12);
    strain = zeros(size(Tet,1),6);
    stress = zeros(size(Tet,1),6);
    vm     = zeros(size(Tet,1),1);
    for i = 1:size(Tet,1)
        B = squeeze(Bs(i,:,:));
        Teti = Tet(i,:);
        ij2p(1:3:end) = 3*Teti-2;
        ij2p(2:3:end) = 3*Teti-1;
        ij2p(3:3:end) = 3*Teti;
        strain(i,:) = B*u(ij2p);
        stress(i,:) = C*strain(i,:)';
        vm(i) = vonmises(stress(i,:));
    end
end
