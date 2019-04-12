% Restores a model
%
function [meshname,V,T,F,b,load] = model_load(inputmodel)
    model = matfile(inputmodel);
    meshname = model.meshname;
    V = model.V;
    T = model.T;
    F = model.F;
    b = model.b;
    load = model.load;
    fprintf('Model (%s)\n#V: %d\n#T: %d\n#b: %d\n', meshname, size(V,1), size(T,1), size(b,1));
end