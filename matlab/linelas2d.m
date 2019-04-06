% 2d linear elasticity
%   reference: https://www.cambro.umu.se/access/content/group/69cbe5bd-67a8-41d8-9bc0-be47f2291b61/B5.pdf
%   requires: PDe toolbox
clear all;

model = createpde(1);
geometryFromEdges(model,@squareg); % dim [-1,1] x [-1,1]
mesh = generateMesh(model, 'GeometricOrder', 'linear');

% https://www.mathworks.com/help/pde/ug/pde.femesh-properties.html
p = mesh.Nodes;
t = mesh.Elements;

[K,F]=assemble(p,[],t);

% enforce dirichlet boundary by adding large number to `K_ii`
for i = 1:size(p,2)
    x = p(1,i); y = p(2,i);
    if x < 0.001-1
        K(2*i-1,2*i-1) = 1.e+6;
        K(2*i,2*i) = 1.e+6;
    end
end


u = K \ F;
deformation.ux = u(1:2:end);
deformation.uy = u(2:2:end);

pdeplot(model, ...
    'Deformation', deformation, ...
    'XYData',  deformation.uy, ...
    'ColorMap','jet')


%   t       3x#F        triangle faces
%   p       2x#V        triangle vertex positions
function [K, F] = assemble(p,e,t)
    ndof = 2*size(p,2);
    K = sparse(ndof,ndof);
    F = zeros(ndof,1);
    dofs = zeros(6,1);
    E = 1; nu = 0.3;
    lambda = E*nu/((1+nu)*(1-2*nu)); mu = E/(2*(1+nu));

    for i = 1:size(t,2)
        nodes = t(1:3,i);
        x = p(1, nodes);
        y = p(2, nodes);
        % indices to put element stiffness to global stiffness matrix
        dofs(1:2:end) = 2*nodes-1;
        dofs(2:2:end) = 2*nodes;
        Ke = stiffness(x,y,mu,lambda);
        Fe = load(x,y);
        K(dofs,dofs) = K(dofs,dofs) + Ke;
        F(dofs) = F(dofs) + Fe;
    end
end


% element stiffness for 2d triangle 
%       x       3x1     x coordinate of 3 node of a triangle
%       y       3x1     y coordinate of 3 node of a triangle
%       mu,lambda   material constants
function Ke = stiffness(x,y,mu,lambda)
    area = polyarea(x,y);
    b = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    c = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    D = mu*[2 0 0; 0 2 0; 0 0 1] + lambda*[1 1 0; 1 1 0; 0 0 0];
    Be = [
        b(1) 0 b(2) 0 b(3) 0;
        0 c(1) 0 c(2) 0 c(3);
        c(1) b(1) c(2) b(2) c(3) b(3)
    ]/2/area;
    Ke = Be'*D*Be*area;
end


function Fe = load(x,y)
    area = polyarea(x,y);
    f = force(mean(x), mean(y));
    Fe = (f(1)*[1 0 1 0 1 0]' + f(2)*[0 1 0 1 0 1]')*area/3;
end


% force at position x,y
function f = force(x,y)
    f = [0 -1];
end