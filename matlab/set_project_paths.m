function set_project_paths()


    % path_to_casadi = '/Users/zhangjiayi/Documents/github/casadi';
    % casadi_subdirs = split(genpath(path_to_casadi),':');
    % addpath(strjoin(gp_subdirs(~contains(casadi_subdirs,'.git')),':'));
    
    % path_to_gptoolbox = '/Users/zhangjiayi/Documents/github/gptoolbox';
    path_to_gptoolbox = '~/github/gptoolbox';
    gp_subdirs = split(genpath(path_to_gptoolbox),':');
    addpath(strjoin(gp_subdirs(~contains(gp_subdirs,'.git')),':'));
    
    % sub-directories
    addpath('.');
    addpath('tests');
    addpath('matrix');
    addpath('mesh');
    addpath('linelas');
    addpath('solvers');
    addpath('utils');
    addpath('scripts');

end