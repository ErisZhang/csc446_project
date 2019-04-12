clear all;
format shortG;
warning('off','all');
set_project_paths;

%% Gen ...

% outputdir = '../../data/outputs/';
% inputmesh = '../data/archbridge_tiny.mesh';
% pathtomodel = '../data/models/archbridge_tiny';
% model_generate(inputmesh,pathtomodel,0.05);
% compute_and_save(outputdir, pathtomodel);






options.face_vertex_color = VM;
plot_mesh(V,F,options);
colormap('jet');

