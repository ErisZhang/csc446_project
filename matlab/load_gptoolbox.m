function load_gptoolbox()
    path_to_gptoolbox = '/Users/markwang/github/gptoolbox/';
    gp_subdirs = split(genpath(path_to_gptoolbox),':');
    addpath(strjoin(gp_subdirs(~contains(gp_subdirs,'.git')),':'));
end