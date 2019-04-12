clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');

b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise


[W,BC,DV,Q,r] = voxelize(V,F,20);


[P,dof,b] = index_ijk_to_p(W);


% d = render_in_cage(V,F,DV,Q,'ColorIntersections',false);


% [Vi,IM] = remove_unreferenced(V,Q);
% Q = IM(Q);
% trisurf(Q,Vi(:,1),Vi(:,2),Vi(:,3));
% axis equal;

load = [0; -9.8; 0];

[U_interp,strain_interp,stress_interp,VM_interp] = linelas3d_hexahedron(W,load,r,DV,V,T);


for i = 1:size(VM_interp,1)
    if VM_interp(i) > 0.75*max(VM_interp)
        VM_interp(i) = 0.75*max(VM_interp);
    end
end

options.face_vertex_color = VM_interp;
plot_mesh(V,F,options);
colormap('jet');

% 
% U_interp = reshape(U_interp.',[],1);
% U = reshape(U.',[],1);
% 
% IJVV = matdiff(U_interp,U,1e-5);
% 
% size(IJVV);
% 
% IJVV(1:10,:)
% 


% FF = zeros(2*size(Q,1),3);
% for i = 1:size(Q,1)
%     idx1 = Q(i,1:3);
%     idx2 = [Q(i,3) Q(i,4) Q(i,1)];
%     FF(2*i-1,:) = idx1;
%     FF(2*i,:) = idx2;
% end
% 
% 
% plot_mesh(DV,FF);
% colormap('jet');

% 
% 
% [P,dof,b] = index_ijk_to_p(W);


% 
% [U_interp] = interpolate_hex_to_tet(P,U,BC,r,V);
