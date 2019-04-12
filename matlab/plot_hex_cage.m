clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/archbridge.mesh');

imgdir = "./data/";
soropt_imageoutput = strcat(imgdir,'hex/');
status = mkdir(soropt_imageoutput);
meshname = "archbridge_cage";

M = [10,20,30,40,60];

for t = 1:5 % 1:20:600
    t;
    clf;
    [W,BC,DV,Q,r] = voxelize(V,F,M(t));
    d = render_in_cage(V,F,DV,Q,'ColorIntersections',false);
    saveas(gcf,strcat(soropt_imageoutput,meshname,'_',num2str(t),'.png'));
end