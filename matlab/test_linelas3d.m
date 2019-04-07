clear all;
warning('off','all');
set_project_paths();

% archbridge.mesh       31,164 tets
[V,T,F] = readMESH('../data/archbridge.mesh');

% U = linelas3d_tetrahedron(V,T);
% save('../data/archbridge_U.mat','U');
bridge = matfile('../data/archbridge_U.mat');
U = bridge.U;






options.face_vertex_color = U(:,1);
plot_mesh(V,F,options);



