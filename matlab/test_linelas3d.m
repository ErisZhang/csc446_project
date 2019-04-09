clear all;
format shortG;
warning('off','all');
set_project_paths();

saved_mats = '../data/archbridge_withbc.mat';
[V,T,F] = readMESH('../data/archbridge.mesh');
b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];

[U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load);
save(saved_mats,'U', 'K', 'f', 'strain', 'stress', 'VM');
% 
bridge = matfile(saved_mats);
U = bridge.U;
K = bridge.K;
f = bridge.f;
strain = bridge.strain;
stress = bridge.stress;
VM = bridge.VM;
u = K\f;

% some hand-tweaking to get color right
VMC = VM;
cap = 0.75*max(VM);
for i = 1:size(VM,1)
    if VMC(i) > cap;
        VMC(i) = cap;
    end
end

options.face_vertex_color = VMC;
plot_mesh(V,F,options);
colormap('jet');
