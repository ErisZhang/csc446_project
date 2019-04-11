% Compute deformation of a 3d solid based on 
%       a linear elasticity model for homogeneous isotropic material
%       where we have hexahedron elements
%
%   Reference: http://what-when-how.com/the-finite-element-method/fem-for-3d-solids-finite-element-method-part-2/
%
% Inputs:
%       V   #Vx3
%       Hex   #Hexx6
%       b   #bx1    boundary index into V
%       len length of the cube
% Outputs:
%       U           #Vx3 list of vertex displacements
%       strain      #Vx6 stress field
%       stress      #Vx6 strain field     
function [U_interp,strain,stress,VM] = linelas3d_hexahedron(W,load,r,BC,V,Tet)
% function [u] = linelas3d_hexahedron(W,load,r,BC,V,Tet)

    assert(size(size(W),2) == 3,'Only 3D meshes are supported');

    young = 1.45e5;
    mu = 0.45;

    len = r(1);
    [P,dof,b] = index_ijk_to_p(W);
    b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1);

    C = [
        (1-mu) mu mu 0 0 0
        mu (1-mu) mu 0 0 0
        mu mu (1-mu) 0 0 0
        0 0 0 0.5*(1-2*mu) 0 0
        0 0 0 0 0.5*(1-2*mu) 0
        0 0 0 0 0 0.5*(1-2*mu)
    ]*(young/(1+mu)/(1-2*mu));


    dim = size(size(W),2);
    dof = dof*dim;

    K = sparse(dof,dof);
    f = zeros(dof,1);

    % vertex postions
    Ve = [-1 -1 -1 % 1/8*(1-i)*(1-j)*(1-k)
           1 -1 -1 % 1/8*(1+i)*(1-j)*(1-k)
           1  1 -1 % 1/8*(1+i)*(1+j)*(1-k)
          -1  1 -1 % 1/8*(1-i)*(1+j)*(1-k)
          -1 -1  1 % 1/8*(1-i)*(1-j)*(1+k)
           1 -1  1 % 1/8*(1+i)*(1-j)*(1+k)
           1  1  1 % 1/8*(1+i)*(1+j)*(1+k)
          -1  1  1];% 1/8*(1-i)*(1+j)*(1+k)

    % gaussian quadrature
    Ke = zeros(24,24);
    weights = [(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 ...
        (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
    points = [-sqrt(5+2*sqrt(10/7))/3 -sqrt(5-2*sqrt(10/7))/3 0 ...
        sqrt(5-2*sqrt(10/7))/3 sqrt(5+2*sqrt(10/7))/3];
    
    N = 5; % number of points for quadrature
    for i = 1:N
        for j = 1:N
            for k = 1:N
                w = weights(i)*weights(j)*weights(k);
                Ke = Ke + w.*f_integrand(C,Ve,points(i),points(j),points(k));
            end
        end
    end

    hex_vol = 2*2*2;
    fe = repmat(load,8,1)*hex_vol/8;

    ij2p = zeros(24,1);
    % assembly procedure
    for i = 1:size(W,1)
        for j = 1:size(W,2)
            for k = 1:size(W,3)
                if W(i,j,k) == 1
                    Hexi = get_hex_index(P,i,j,k);
                    ij2p(1:3:end) = 3*Hexi-2;
                    ij2p(2:3:end) = 3*Hexi-1;
                    ij2p(3:3:end) = 3*Hexi;

                    K(ij2p,ij2p) = K(ij2p,ij2p) + Ke;
                    f(ij2p) = f(ij2p) + fe;
                end
            end
        end
    end

    [K,f] = dirichlet_zero_boundary(K,f,b);

    u = K \ f;
    
    u = u.*len.*len./4;



    U = zeros(dof/dim,dim);
    U(:,1) = u(1:3:end);
    U(:,2) = u(2:3:end);
    U(:,3) = u(3:3:end);

    [U_interp] = interpolate_hex_to_tet(P,U,BC,r,V);
    u_interp = zeros(size(V,2)*size(V,1),1);
    u_interp(1:3:end) = U_interp(:,1);
    u_interp(2:3:end) = U_interp(:,2);
    u_interp(3:3:end) = U_interp(:,3);

    B = zeros(6,12);
    Bs = zeros(size(Tet,1),6,12);

    for i = 1:size(Tet,1)

        Teti = Tet(i,:);

        % vertex positions
        v1 = V(Teti(1),:);
        v2 = V(Teti(2),:);
        v3 = V(Teti(3),:);
        v4 = V(Teti(4),:);

        % barycentric -> cartesian
        T = [
            1 1 1 1
            v1(1) v2(1) v3(1) v4(1)
            v1(2) v2(2) v3(2) v4(2)
            v1(3) v2(3) v3(3) v4(3)
        ];

        % volume of a tetrahedron
        tet_vol = (1/6)*det(T);

        % cartesian -> barycentric
        Tinv = inv(T)*6*tet_vol;

        % B matrix
        %       where \epsilon = B * u^e
        %       B   6x12
        %       u   12x1
        B(1,1:3:12) = Tinv(:,2);
        B(2,2:3:12) = Tinv(:,3);
        B(3,3:3:12) = Tinv(:,4);
        B(4,:) = circshift(B(1,:),1) + circshift(B(2,:),-1);
        B(5,:) = circshift(B(3,:),-1) + circshift(B(2,:),1);
        B(6,:) = circshift(B(3,:),-2) + circshift(B(1,:),2);
        B = B/(6*tet_vol);
        Bs(i,:,:)=B;
    end

    [strain,stress,vm]=per_element_fields(Tet,C,Bs,u_interp);

    N = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        N(Tet(i,:)) = N(Tet(i,:)) + 1;
    end

    VM = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        VM(Tet(i,:),1) = VM(Tet(i,:),1) + vm(i)./N(Tet(i,:),1);
    end


end


function Hexi = get_hex_index(P,i,j,k)
    Hexi = zeros(8,1);
    for a = 0:1
        for b = 0:1
            for c = 0:1
                idx = 4*a+2*b+c+1;
                Hexi(idx) = P(i+a,j+b,k+c);
            end
        end
    end
end


function value = f_integrand(C,Ve,i,j,k)

    J = [-1/8*(1-j)*(1-k)  1/8*(1-j)*(1-k)  1/8*(1+j)*(1-k) -1/8*(1+j)*(1-k) -1/8*(1-j)*(1+k)  1/8*(1-j)*(1+k) 1/8*(1+j)*(1+k) -1/8*(1+j)*(1+k)
         -1/8*(1-i)*(1-k) -1/8*(1+i)*(1-k)  1/8*(1+i)*(1-k)  1/8*(1-i)*(1-k) -1/8*(1-i)*(1+k) -1/8*(1+i)*(1+k) 1/8*(1+i)*(1+k)  1/8*(1-i)*(1+k)
         -1/8*(1-i)*(1-j) -1/8*(1+i)*(1-j) -1/8*(1+i)*(1+j) -1/8*(1-i)*(1+j)  1/8*(1-i)*(1-j)  1/8*(1+i)*(1-j) 1/8*(1+i)*(1+j)  1/8*(1-i)*(1+j)]*Ve;  
    
    B = [];
    for m = 1:8
        coeff = Ve(m,:);
        
        Ni_xyz = inv(J)*[1/8*coeff(1)*(1+coeff(2)*j)*(1+coeff(3)*k)
                         1/8*coeff(2)*(1+coeff(1)*i)*(1+coeff(3)*k)     
                         1/8*coeff(3)*(1+coeff(1)*i)*(1+coeff(2)*j)];
        
        % assemble matrix B
        B = [B [Ni_xyz(1)         0         0
                        0 Ni_xyz(2)         0
                        0         0 Ni_xyz(3)
                        0 Ni_xyz(3) Ni_xyz(2)
                Ni_xyz(3)         0 Ni_xyz(1)
                Ni_xyz(2) Ni_xyz(1)         0]];

    end

    value = B'*C*B*det(J);

end