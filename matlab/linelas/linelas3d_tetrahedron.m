% Compute deformation of a 3d solid based on 
%       a linear elasticity model for homogeneous isotropic material
%       where we have tetrahedron elements
%
%       reference: http://automatica.dei.unipd.it/public/Schenato/PSC/2010_2011/gruppo4-Building_termo_identification/IdentificazioneTermodinamica20062007/eggpsc/Simulatore%20termico/AFEM.Ch15.pdf
% Inputs:
%       V   #Vx2
%       T   #Tetx3
% Outputs:
%       U   #Vx2 list of vertex displacements
function U = linelas3d_tetrahedron(V,Tet)
    assert(size(V,2) == 3,'Only 3D meshes are supported');

    a = 1e6;    % young modulus
    b = 0.45;   % poisson ratio

    C = [
        (1-b) b b 0 0 0
        b (1-b) b 0 0 0
        b b (1-b) 0 0 0
        0 0 0 0.5*(1-2*b) 0 0
        0 0 0 0 0.5*(1-2*b) 0
        0 0 0 0 0 0.5*(1-2*b)
    ]*(a/(1+b)/(1-2*b));


    load = [-9.8; 0; 0];

    dim = size(V,2);
    dof = size(V,1)*dim;

    K = sparse(dof,dof);
    F = zeros(dof,1);
    
    ij2p = zeros(12,1);
    B = zeros(6,12);
    Ke = zeros(12,12);
    Fe = zeros(12,1);
    T = zeros(4,4);
    Tinv = zeros(4,4);

    tic;
    for i = 1:size(Tet,1)

        if mod(i,1000) == 0
            [i size(Tet,1)]
        end

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
        Tinv = inv(T);

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

        % element stiffness
        %       12x12
        Ke = tet_vol*B'*C*B;

        % body force 
        %       12x1
        Fe = repmat(load,4,1)*tet_vol/4;

        % assembly

        % local -> global indexing 
        ij2p(1:3:end) = 3*Teti-2;
        ij2p(2:3:end) = 3*Teti-1;
        ij2p(3:3:end) = 3*Teti;

        K(ij2p,ij2p) = K(ij2p,ij2p) + Ke;
        F(ij2p) = F(ij2p) + Fe;
    end
    toc;    % 68s

    % % tested once \checkmark
    % assert(issymmetric_approx(K, 1e-3) == true);
    % assert(ispd(K) == true);

    tic;
    u = K \ F;
    toc;    % 1.45s

    U = zeros(size(V));
    U(:,1) = u(1:3:end);
    U(:,2) = u(2:3:end);
    U(:,3) = u(3:3:end);
end