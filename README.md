


```
# clone https://github.com/alecjacobson/gptoolbox
gp_subdirs = split(genpath('/path/to/gptoolbox'),':');
addpath(strjoin(gp_subdirs(~contains(gp_subdirs,'.git')),':'));


```