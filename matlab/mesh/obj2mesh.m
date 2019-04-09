% need to instal cmd tool `tetgen`
%   and compile mex https://github.com/alecjacobson/gptoolbox/tree/master/mex
%   
%   If ratio \neq 1, use the ratio for decimate ...
%   
function obj2mesh(objfile,meshfile,ratio)
    [V,F,~,~,~,~] = readOBJ(objfile);
    if ratio ~= 1
        [V,F] = decimate_libigl(V,F,ratio);
    end
    [V,T,F] = tetgen(V,F,'Verbose',true,'Flags','-q100');
    writeMESH(meshfile,V,T,F,[]);
end