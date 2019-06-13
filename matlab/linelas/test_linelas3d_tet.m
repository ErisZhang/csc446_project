clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/bb-bunny.mesh');


b = boundary_vertices(V,2,0.005);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];

[U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load);

options.face_vertex_color = log(VM);
plot_mesh(V,F,options);
colormap('jet');