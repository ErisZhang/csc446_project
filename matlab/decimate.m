clear all;
format shortG;
warning('off','all');
set_project_paths();


obj2mesh("/Users/zhangjiayi/Documents/csc446_project/data/armadillo.obj",...
    "/Users/zhangjiayi/Documents/csc446_project/data/armadillo.mesh",0.1);