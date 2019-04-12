clear all;
format shortG;
set_project_paths();

[V,T,F] = readMESH('../data/archbridge_tiny.mesh');

load = [0; -9.8; 0];

max_iters = 300;
tol = 1e-20;
omega =  1.5;
saveon = (1:5:300)';

[W,BC,DV,Q,r] = voxelize(V,F,40);

[U_interp,strain_interp,stress_interp,VM,P,C,data] = linelas3d_hexahedron(W,load,r,DV,V,T, ...
    'LinearSolver', @(A,b) sor(A,b,zeros(size(b)),max_iters,tol,'Omega',omega,'SaveOn',saveon));

% L2 column (relative residual) norm over data.rks
% plot(1:size(data.rks,2), log(vecnorm(data.rks) ./ norm(b)));

for t = 1:size(data.xks,2)
    
    xk = data.xks(:,t);
    
    [U_interp,strain,stress,VMt] = compute_fields_hex(P,DV,V,T,C,xk,r);
    
    if t == size(data.xks,2)
        options.face_vertex_color = VMt;
        plot_mesh(V,F,options);
        colormap('jet');
        mean((VMt-VM)./VMt)
    end
    pause(0.1)
end


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