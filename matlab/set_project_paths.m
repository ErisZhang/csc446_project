function set_project_paths()

    % gptoolbox 
    path_to_gptoolbox = '/Users/markwang/github/gptoolbox/';
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