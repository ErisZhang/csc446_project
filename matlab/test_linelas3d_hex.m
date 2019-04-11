clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');
b = boundary_vertices(V,2,0.02);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise


[W,BC,DV,Q,r] = voxelize(V,F,30);


[V,Q] = voxel_surface(W,'Centers',DV);


% [Vi,IM] = remove_unreferenced(V,Q);
% Q = IM(Q);
% trisurf(Q,Vi(:,1),Vi(:,2),Vi(:,3));
% axis equal;

% load = [0; -9.8; 0];
% 
% [U_interp,strain_interp,stress_interp,VM_interp] = linelas3d_hexahedron(W,load,r,BC,V,T);

% [U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load);
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



FF = zeros(2*size(Q,1),3);
for i = 1:size(Q,1)
    idx1 = Q(i,1:3);
    idx2 = [Q(i,3) Q(i,4) Q(i,1)];
    FF(2*i-1,:) = idx1;
    FF(2*i,:) = idx2;
end


plot_mesh(V,FF);
colormap('jet');

% 
% 
% [P,dof,b] = index_ijk_to_p(W);


% 
% [U_interp] = interpolate_hex_to_tet(P,U,BC,r,V);
% 
% [DV,I] = remove_unreferenced(DV,Q);
% Q = I(Q);
% d = render_in_cage(V,F,DV,Q,'ColorIntersections',false);
