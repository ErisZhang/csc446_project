% adds boundary dirichlet zero boundary condition to `K`, `f`
%       where `K` is sparse and
%             'b' is index to `K`/`f` for enforcing zero boundary
%                   can specify boundary for arbitrary vertex, 
%                   in arbitrary direction {x,y,z}
%
%       note: for 3d problem where `K`/`f` is node-wise indexed,
%               
%           b = {vertex indices};
%           b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1)
%   
function [Kz,fz] = dirichlet_zero_boundary(K,f,b)
    Kz = K;
    fz = f;
    
    [I,J,~]=find(K);
    for i=1:size(I,1)
        if any(I(i) == b) || any(J(i) == b)
            Kz(I(i),J(i)) = 0;
        end
    end
    for i=1:size(b)
        Kz(b(i),b(i)) = 1;
        fz(b(i)) = 0;
    end
end
    