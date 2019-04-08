% need to instal cmd tool `tetgen`
function obj2mesh(objfile,meshfile)
    [V,F,~,~,~,~] = readOBJ(objfile);
    [V,T,F] = tetgen(V,F,'Verbose',true,'Flags','-q100');
    writeMESH(meshfile,V,T,F,[]);
end