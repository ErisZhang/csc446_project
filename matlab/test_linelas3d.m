clear all;
format shortG;
warning('off','all');
set_project_paths();

saved_mats = '../data/archbridge_withbc.mat';
[V,T,F] = readMESH('../data/archbridge_tiny.mesh');
b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];

[U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load)
% [U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load, ...
%     'LinearSolver', @(A,b) sor(A,b,zeros(size(b)),1e3,1e-11,sor_omegaopt(A)));
% save(saved_mats,'U', 'K', 'f', 'strain', 'stress', 'VM');

bridge = matfile(saved_mats);
U = bridge.U;
K = bridge.K;
f = bridge.f;
strain = bridge.strain;
stress = bridge.stress;
VM = bridge.VM;
cond(K)

% compare linear solvers ..
% u = K\f;
% uhat = reshape(U.',[],1);
% IJVV = matdiff(u,uhat,1e-5);
% size(IJVV)
% IJVV(1:10,:)


options.face_vertex_color = VM;
plot_mesh(V,F,options);
colormap('jet');

