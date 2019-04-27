clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');


[W,BC,DV,Q,r] = voxelize(V,F,40);


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
