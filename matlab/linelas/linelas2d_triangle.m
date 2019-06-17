function [U,strain,stress,vm,VM] = linelas2d_triangle(V,F,b,load)
    young = 1.45e5;
    mu = 0.45;

    D = [1 mu 0
        mu  1 0
         0  0  0.5*(1-mu)]*young/(1-mu^2);

    dim = size(V,2); % vertices live in 2d
    dof = size(V,1)*dim; % number of vertices

    % global stiffness matrix
    K = sparse(dof,dof);
    f = zeros(dof,1);

    % element stiffness matrix
    B = zeros(3,6);
    Ke = zeros(6,6);
    fe = zeros(6,1);
    ij2p = zeros(6,1);

    % store B and Vo for each element
    Bs = zeros(size(F,1),3,6);
    Vs = zeros(size(F,1),1);


    for i = 1:size(F,1)

        ele_i = F(i,:);

        % vertex positions
        v1 = V(ele_i(1),:);
        v2 = V(ele_i(2),:);
        v3 = V(ele_i(3),:);

        C = [1 v1
             1 v2
             1 v3];

        IC = inv(C);

        tri_area = det(C)/2;

        for j=1:3
            B(1,2*j-1)=IC(2,j);
            B(1,2*j)=0;
            B(2,2*j-1)=0;
            B(2,2*j)=IC(3,j);
            B(3,2*j-1)=IC(3,j);
            B(3,2*j)=IC(2,j);
        end

        Ke = B'*D*B*tri_area;

        Bs(i,:,:) = B;

        fe = repmat(load,3,1)*tri_area/3;

        % local -> global indexing 
        ij2p(1:2:end) = 2*ele_i-1;
        ij2p(2:2:end) = 2*ele_i-0;

        K(ij2p,ij2p) = K(ij2p,ij2p) + Ke;
        f(ij2p) = f(ij2p) + fe;

    end

    [K,f] = dirichlet_zero_boundary(K,f,b);

    u = K\f;
    
    [U,strain,stress,vm,VM] = compute_fields_2d(V,F,D,Bs,u);

end