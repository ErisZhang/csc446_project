clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/archbridge.mesh');

imgdir = "./data/";
soropt_imageoutput = strcat(imgdir,'hex/');
status = mkdir(soropt_imageoutput);
meshname = "archbridge";

M = [10,20,30,40,60];

for t = 1:5 % 1:20:600
    t;
    clf;
    [W,BC,DV,Q,r] = voxelize(V,F,M(t));
    load = [0; -9.8; 0];
    [U_interp,strain_interp,stress_interp,VM_interp] = linelas3d_hexahedron(W,load,r,DV,V,T);
    for i = 1:size(VM_interp,1)
        if VM_interp(i) > 0.78*max(VM_interp)
            VM_interp(i) = 0.78*max(VM_interp);
        end
    end
    load = [0; -9.8; 0];
    options.face_vertex_color = VM_interp;
    plot_mesh(V,F,options);
    colormap('jet');
    view(90,-50);
    camroll(90);
    shading interp;
    axis tight;
    saveas(gcf,strcat(soropt_imageoutput,meshname,'_',num2str(t),'.png'));
end