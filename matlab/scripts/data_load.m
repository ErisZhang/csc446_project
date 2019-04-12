% load precomputed data under `data/output`
function [pathtomodel,solvername,omega,max_iters,tol,saveon,U,K,f,strain,stress,VM,Bs,C,data] = data_load(matname)
    d = matfile(matname);
    pathtomodel     =       d.pathtomodel;
    solvername      =       d.solvername;
    omega       =       d.omega;
    max_iters       =       d.max_iters;
    tol     =       d.tol;
    saveon      =       d.saveon;
    U       =       d.U;
    K       =       d.K;
    f       =       d.f;
    strain      =       d.strain;
    stress      =       d.stress;
    VM      =       d.VM;
    Bs      =       d.Bs;
    C       =       d.C;
    data        =       d.data;
    fprintf('Data (%s)\nsolvername: %s\nomega: %f\nsize(K): %d x %d \n',matname,solvername,omega,size(K,1),size(K,2));
end
