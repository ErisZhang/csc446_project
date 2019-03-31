clear all;
warning('off','all');
load_gptoolbox();

[V,T,F] = readMESH('/Users/markwang/github/fast_support_reduction/data/archbridge.mesh');
C = zeros(size(V));
C(:,1:1) = 1;
options.face_vertex_color = jet(size(V,1));
plot_mesh(V,F,options);