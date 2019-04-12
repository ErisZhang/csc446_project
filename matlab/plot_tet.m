clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/archbridge.mesh');

imgdir = "./data/";
soropt_imageoutput = strcat(imgdir,'tet/');
status = mkdir(soropt_imageoutput);
meshname = "archbridge";

t = 1;

b = boundary_vertices(V,2,0.005);
b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
load = [0; -9.8; 0];

[U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b,load);

for i = 1:size(VM,1)
    if VM(i) > 0.85*max(VM)
        VM(i) = 0.85*max(VM);
    end
end

options.face_vertex_color = VM;

plot_mesh(V,F,options);
colormap('jet');
view(90,-50);
camroll(90);
shading interp;
axis tight;

saveas(gcf,strcat(soropt_imageoutput,meshname,'_',num2str(t),'.png'));