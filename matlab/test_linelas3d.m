clear all;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/archbridge.mesh');
C = zeros(size(V));
C(:,1:1) = 1;
C = jet(size(V,1))
options.face_vertex_color = C;
plot_mesh(V,F,options);

