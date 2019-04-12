clear all;
format shortG;
warning('off','all');
set_project_paths;

imgdir = '../data/imgs/';
status = mkdir(imgdir);
outputdir = '../data/outputs/';
inputmesh = '../data/archbridge_tiny.mesh';
pathtomodel = '../data/models/archbridge_tiny';
[meshname,V,T,F,bnds,load] = model_load(pathtomodel);



%%%%%%%%%%%%%%%%%%%%%
% Load SOR into memory
%%%%%%%%%%%%%%%%%%%%%

solvername = 'sor';
omegas = 1:0.2:1.8;
inmem = {};
for it = 1:size(omegas,2)
    omega = omegas(it);
    matname = strcat(outputdir,meshname,'_',solvername,num2str(omega),'.mat');
    [pathtomodel,solvername,omega,max_iters,tol,saveon,U,K,f,strain,stress,VM,Bs,C,data] = data_load(matname);

    %
    inmem{it}.solvername = solvername;
    inmem{it}.saveon = saveon;
    inmem{it}.omega = omega;
    inmem{it}.max_iters = max_iters;
    inmem{it}.tol = tol;
    inmem{it}.U = U;
    inmem{it}.K = K;
    inmem{it}.f = f;
    inmem{it}.VM = VM;
    inmem{it}.Bs = Bs;
    inmem{it}.C = C;
    inmem{it}.data = data;
    %
end


%%%%%%%%%%%%%%%%%%%%%
% SOR relative residual vs. iterations
%%%%%%%%%%%%%%%%%%%%%

hold on;
whichks = 1:1:400;
names = {};
for it = 1:size(omegas,2)
    %
    saveon = inmem{it}.saveon;
    data = inmem{it}.data;
    omega = inmem{it}.omega;
    f = inmem{it}.f;
    %

    ks = saveon(whichks);
    % relres = log(vecnorm(data.rks) ./ norm(f));
    relres = vecnorm(data.rks) ./ norm(f);
    relres = relres(whichks);
    plot(ks, relres, 'LineWidth', 1.5);
    names{it} = strcat('SOR (ω = ',num2str(omega),')');
end
legend(names);
grid on;
xlabel('Number of Iterations')
ylabel('Relative Residual')
saveas(gcf,strcat(imgdir,meshname,'_',solvername,'_relres_vs_iters.png'));
hold off;



%%%%%%%%%%%%%%%%%%%%%
% SOR log relative residual vs. iterations
%%%%%%%%%%%%%%%%%%%%%

hold on;
whichks = 1:1:1000;
names = {};
for it = 1:size(omegas,2)
    %
    saveon = inmem{it}.saveon;
    data = inmem{it}.data;
    omega = inmem{it}.omega;
    f = inmem{it}.f;
    %

    ks = saveon(whichks);
    relres = log(vecnorm(data.rks) ./ norm(f));
    relres = relres(whichks);
    plot(ks, relres, 'LineWidth', 1.5);
    names{it} = strcat('SOR (ω = ',num2str(omega),')');
end
legend(names);
grid on;
xlabel('Number of Iterations')
ylabel('Log Relative Residual')
saveas(gcf,strcat(imgdir,meshname,'_',solvername,'_logrelres_vs_iters.png'));
hold off;



%%% hack

sormem = inmem;

%%%%%%%%%%%%%%%%%%%%%
% Load Jacobi into memory 
%%%%%%%%%%%%%%%%%%%%%


solvername = 'jacobi';
omegas = 0:0.08:0.4;
omegas = [omegas 0.6 0.8 1]
jacobimem = {};
for it = 1:size(omegas,2)
    omega = omegas(it);
    matname = strcat(outputdir,meshname,'_',solvername,num2str(omega),'.mat');
    [pathtomodel,solvername,omega,max_iters,tol,saveon,U,K,f,strain,stress,VM,Bs,C,data] = data_load(matname);

    %
    jacobimem{it}.solvername = solvername;
    jacobimem{it}.saveon = saveon;
    jacobimem{it}.omega = omega;
    jacobimem{it}.max_iters = max_iters;
    jacobimem{it}.tol = tol;
    jacobimem{it}.U = U;
    jacobimem{it}.K = K;
    jacobimem{it}.f = f;
    jacobimem{it}.VM = VM;
    jacobimem{it}.Bs = Bs;
    jacobimem{it}.C = C;
    jacobimem{it}.data = data;
    %
end


%%%%%%%%%%%%%%%%%%%%%
% Jacobi relative residual vs. iterations
%%%%%%%%%%%%%%%%%%%%%

hold on;
whichks = 1:1:1000;
names = {};
for it = 1:size(omegas,2)
    %
    saveon = jacobimem{it}.saveon;
    data = jacobimem{it}.data;
    omega = jacobimem{it}.omega;
    f = jacobimem{it}.f;
    %

    ks = saveon(whichks);
    % relres = log(vecnorm(data.rks) ./ norm(f));
    relres = vecnorm(data.rks) ./ norm(f);
    relres = relres(whichks);
    plot(ks, relres, 'LineWidth', 1.5);
    names{it} = strcat('Jacobi (ω = ',num2str(omega),')');
end
legend(names);
grid on;
ylim([0 1])
xlabel('Number of Iterations')
ylabel('Relative Residual')
saveas(gcf,strcat(imgdir,meshname,'_',solvername,'_relres_vs_iters_full.png'));
hold off;




%%%%%%%%%%%%%%%%%%%%%
% Jacobi Log relative residual vs. iterations
%%%%%%%%%%%%%%%%%%%%%

hold on;
whichks = 1:1:1000;
names = {};
for it = 1:size(omegas,2)
    %
    saveon = jacobimem{it}.saveon;
    data = jacobimem{it}.data;
    omega = jacobimem{it}.omega;
    f = jacobimem{it}.f;
    %

    ks = saveon(whichks);
    relres = log(vecnorm(data.rks) ./ norm(f));
    relres = relres(whichks);
    plot(ks, relres, 'LineWidth', 1.5);
    names{it} = strcat('Jacobi (ω = ',num2str(omega),')');
end
legend(names);
grid on;
ylim([-1 0.1])
xlabel('Number of Iterations')
ylabel('Log Relative Residual')
saveas(gcf,strcat(imgdir,meshname,'_',solvername,'_logrelres_vs_iters_full.png'));
hold off;



%%%%%%%%%%%%%%%%%%%%%
% Jacobi: w = 0.4
% SOR:    w = 1.6 (1.8 long run)
%%%%%%%%%%%%%%%%%%%%%

wjac = 0.4;
wsor = 1.6;

jacopt = jacobimem{6};
soropt = sormem{4};


%%%%%%%%%%%%%%%%%%%%%
% Load goldstandard data using backslash
%%%%%%%%%%%%%%%%%%%%%



[pathtomodel,solvername,omega,max_iters,tol,saveon,U,K,f,strain,stress,VM,Bs,C,data] = data_load(strcat(outputdir,meshname,'_','backslash.mat'));
bsmem.solvername = solvername;
bsmem.saveon = saveon;
bsmem.omega = omega;
bsmem.max_iters = max_iters;
bsmem.tol = tol;
bsmem.U = U;
bsmem.K = K;
bsmem.f = f;
bsmem.VM = VM;
bsmem.Bs = Bs;
bsmem.C = C;
bsmem.data = data;



%%%%%%%%%%%%%%%%%%%%%
% Given goldstandard with backslash, compute absolute error norm with jacopt and soropt
%%%%%%%%%%%%%%%%%%%%%

hold on;
whichks = 1:1:400;
names = {
    'Jacobi (ω = 0.4)',
    'SOR    (ω = 1.6)'
};

bsmem_referenceu = reshape(bsmem.U',[],1);
ks = bsmem.saveon(whichks);

jacobi_err = vecnorm(bsmem_referenceu - jacopt.data.xks);
jacobi_err = jacobi_err(whichks);
% jacobi_err = log(jacobi_err);
plot(ks,jacobi_err,'LineWidth', 1.5);
sor_err = vecnorm(bsmem_referenceu - soropt.data.xks);
% sor_err = log(sor_err);
sor_err = sor_err(whichks);
plot(ks, sor_err,'LineWidth', 1.5);

legend(names);
grid on;
% ylim([-1 0.1])
xlabel('Iterations');
ylabel('Error Norm');  % assuming / solution is correct
saveas(gcf,strcat(imgdir,meshname,'_errornorm.png'));
hold off;



%%%%%%%%%%%%%%%%%%%%%
% Log scale of above
%%%%%%%%%%%%%%%%%%%%%


hold on;
whichks = 1:1:400;
names = {
    'Jacobi (ω = 0.4)',
    'SOR    (ω = 1.6)'
};

bsmem_referenceu = reshape(bsmem.U',[],1);
ks = bsmem.saveon(whichks);

jacobi_err = vecnorm(bsmem_referenceu - jacopt.data.xks);
jacobi_err = jacobi_err(whichks);
jacobi_err = log(jacobi_err);
plot(ks,jacobi_err,'LineWidth', 1.5);
sor_err = vecnorm(bsmem_referenceu - soropt.data.xks);
sor_err = log(sor_err);
sor_err = sor_err(whichks);
plot(ks, sor_err,'LineWidth', 1.5);

legend(names);
grid on;
% ylim([-1 0.1])
xlabel('Iterations');
ylabel('Log Error Norm');  % assuming / solution is correct
saveas(gcf,strcat(imgdir,meshname,'_logerrornorm.png'));
hold off;

%%%%%%%%%%%%%%%%%%%%%
% Visualization at different iterations
%%%%%%%%%%%%%%%%%%%%%

jacopt_imageoutput = strcat(imgdir,'jacopt_vm/');
status = mkdir(jacopt_imageoutput);
C = jacopt.C;
Bs = jacopt.Bs;
ts = [1:5:100 100:20:1000]
for t = ts
    t
    clf;
    xk = jacopt.data.xks(:,t);
    [~,~,~,~,VM] = compute_fields(V,T,C,Bs,xk);
    options.face_vertex_color = VM;
    plot_mesh(V,F,options);
    colormap('jet');
    view(90,-50);
    camroll(90);
    shading interp;
    axis tight;
    pause(0.1);
    saveas(gcf,strcat(jacopt_imageoutput,meshname,'_',num2str(t),'.png'));
end



%%%%%%%%%%%%%%%%%%%%%
% Visualization at different iterations for SOR
%%%%%%%%%%%%%%%%%%%%%

soropt_imageoutput = strcat(imgdir,'soropt_vm/');
status = mkdir(soropt_imageoutput);
C = soropt.C;
Bs = soropt.Bs;
ts = [1:5:100 100:20:1000]
for t = ts
    t
    clf;
    xk = soropt.data.xks(:,t);
    [~,~,~,~,VM] = compute_fields(V,T,C,Bs,xk);
    options.face_vertex_color = VM;
    plot_mesh(V,F,options);
    colormap('jet');
    view(90,-50);
    camroll(90);
    shading interp;
    axis tight;
    pause(0.1);
    saveas(gcf,strcat(soropt_imageoutput,meshname,'_',num2str(t),'.png'));
end



%%%%%%%%%%%%%%%%%%%%%
% Visualization for goldstandard 
%%%%%%%%%%%%%%%%%%%%%
bsmem_imageoutput = strcat(imgdir,'goldstandard_vm/');
status = mkdir(bsmem_imageoutput);
C = bsmem.C;
Bs = bsmem.Bs;
bsmem_referenceu = reshape(bsmem.U',[],1);
[~,~,~,~,VM] = compute_fields(V,T,C,Bs,bsmem_referenceu);
options.face_vertex_color = VM;
plot_mesh(V,F,options);
colormap('jet');
view(90,-50);
camroll(90);
shading interp;
axis tight;
saveas(gcf,strcat(bsmem_imageoutput,meshname,'_goldstandard.png'));





% whichks = 1:500
% ks = saveon(whichks);
% relres = log(vecnorm(data.rks) ./ norm(f));
% relres = relres(whichks);




% [V,T,F] = readMESH('../data/archbridge_tiny.mesh');
% bridge = matfile('../data/archbridge_withbc.mat');
% b = boundary_vertices(V,2,0.02);
% b = reshape(3*repmat(b,3,1) - [2 1 0]',[],1); % vertex->node-wise
% load = [0; -9.8; 0];
% U = bridge.U;
% K = bridge.K;
% f = bridge.f;
% VM=  bridge.VM;
% u = K\f;
% condest(K) 

% x0 = zeros(size(f));
% max_iters = 300;
% tol = 1e-20;
% omega =  1.5;
% saveon = (1:5:300)';

saveons = [1;2;3];
linelas3d_tetrahedron(V,T,bnds,load, 'LinearSolver', @(A,b) sor(A,b,zeros(size(b)),10,tol,'Omega',1.8,'SaveOn',saveons));

% for t = 1:size(data.xks,2)
%     xk = data.xks(:,t);
%     [~,~,~,~,VM] = compute_fields(V,T,C,Bs,xk);
%     options.face_vertex_color = VM;
%     plot_mesh(V,F,options);
%     colormap('jet');
%     mean((VM-VMh)./VM)
%     pause(0.2)
% end

% 
% 
% IJVV = matdiff(U,Uh,1e-6);
% size(IJVV)
% IJVV(1:min(10,end),:);

% options.face_vertex_color = VMh;
% plot_mesh(V,F,options);
% colormap('jet');




% options.face_vertex_color = (VM-VMh)./VM;
% plot_mesh(V,F,options);
% colormap(redblue);

% vm stress error

% `jet` by default has 64 colors
% so 1/64=0.0156
% 1/256 = 0.0039