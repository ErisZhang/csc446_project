clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');


[W,BC,DV,Q,r] = voxelize(V,F,40);


load = [0; -9.8; 0];

[VM,P,C,data] = linelas3d_hexahedron(W,load,r,DV,V,T);

ori = DV(1,:);

% for i = 1:size(Q,1)
%     for j = 1:4
%         idx = floor((DV(Q(i,j),:)-ori)./r)+1;
%         idx = index_ijk_to_p(idx);
%         if idx > size(VM)
%             disp("problem here");
%         end
%     end
% end

% 
% options.face_vertex_color = VM;
% plot_mesh(V,F,options);
% colormap('jet');
