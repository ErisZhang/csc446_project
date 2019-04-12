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
function [U_interp,strain,stress,VM,P,C,data] = linelas3d_hexahedron(W,load,r,DV,V,Tet,varargin)
    assert(size(size(W),2) == 3,'Only 3D meshes are supported');

    young = 1.45e5;
    mu = 0.45;
    linearsolver = false;
    saveon = {};
    data = struct('xks',[],'rks',[]);

    % Map of parameter names to variable names
    params_to_variables = containers.Map( ...
        {'Young','Mu','LinearSolver', 'SaveOn'}, ...
        {'young','mu','linearsolver', 'saveon'});
    v = 1;
    while v <= numel(varargin)
        param_name = varargin{v};
        if isKey(params_to_variables,param_name)
            assert(v+1<=numel(varargin));
            v = v+1;
            % Trick: use feval on anonymous function to use assignin to this workspace
            feval(@()assignin('caller',params_to_variables(param_name),varargin{v}));
        else
            error('Unsupported parameter: %s',varargin{v});
        end
        v=v+1;
    end



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

    tic;
    if islogical(linearsolver)
        u = K\f;
    else
        [u,data.xks,data.rks] = linearsolver(K,f);
    end
    toc;
    
    u = u.*len.*len./4;

    [U_interp,strain,stress,VM] = compute_fields_hex(P,DV,V,Tet,C,u,r);


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