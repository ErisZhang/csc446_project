clear all;
format shortG;
warning('off','all');
set_project_paths();

saved_mats = '../data/archbridge_withbc.mat';
[V,T,F] = readMESH('../data/archbridge.mesh');
b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];

% [U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load);
% save(saved_mats,'U', 'K', 'f', 'strain', 'stress', 'VM');
% 
bridge = matfile(saved_mats);
U = bridge.U;
K = bridge.K;
f = bridge.f;
strain = bridge.strain;
stress = bridge.stress;
VM = bridge.VM;

tic;
u = K\f;
toc;

tic;
eps = 1e-11;
uhat = gaussseidels(K,f,zeros(size(f)),1000000000,eps);
toc;

uuhat = [u uhat];
uuhat(1:10,:)

Ku = full(K);
Kudiag = diag(Ku).^-1;
cond(Kudiag*K)



assert(all(abs(u-uhat)./abs(u) < eps*100), 'did not converge ...');

% some hand-tweaking to get color right
VMC = VM;
size(VM)
cap = 0.75*max(VM);
for i = 1:size(VM,1)
    if VMC(i) > cap
        VMC(i) = cap;
    end
end

options.face_vertex_color = VMC;
plot_mesh(V,F,options);
colormap('jet');
