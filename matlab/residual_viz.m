clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');
bridge = matfile('../data/archbridge_withbc.mat');
b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];
U = bridge.U;
K = bridge.K;
f = bridge.f;
VM=  bridge.VM;
u = K\f;
condest(K) 

x0 = zeros(size(f));
max_iters = 3000;
tol = 1e-12;
omega =  1.5;


[Uh,Kh,fh,strainh,stressh,VMh] = linelas3d_tetrahedron(V,T,b,load, ...
    'LinearSolver', @(A,b) sor(A,b,zeros(size(b)),max_iters,tol,'Omega',omega));

IJVV = matdiff(U,Uh,1e-6);
size(IJVV)
IJVV(1:min(10,end),:);

% options.face_vertex_color = VMh;
% plot_mesh(V,F,options);
% colormap('jet');


% options.face_vertex_color = VM;
% plot_mesh(V,F,options);
% colormap('jet');


options.face_vertex_color = (VM-VMh)./VM;
plot_mesh(V,F,options);
colormap(redblue);

% vm stress error
mean((VM-VMh)./VM)

% `jet` by default has 64 colors
% so 1/64=0.0156
% 1/256 = 0.0039