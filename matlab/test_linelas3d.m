clear all;
format shortG;
warning('off','all');
set_project_paths();

% archbridge.mesh       31,164 tets

saved_mats = '../data/archbridge_nobc.mat';
[V,T,F] = readMESH('../data/archbridge.mesh');

% saved_mats = '../data/bunny_computed.mat';
% [V,T,F] = readMESH('../data/bunny.mesh');

% V = V + abs(min(V)); % to 1st quadrant

% boundary vertices
tol = 0.05*abs(max(V(:,2))-min(V(:,2)));
b = []; nbc = 1;
for i = 1:size(V,1)
    if abs(V(i,2)-min(V(:,2))) < tol
        b(nbc) = i;
        nbc = nbc + 1;
    end
end


C = zeros(size(V)) + [0.3,0.3,0.3];
for i = 1:max(size(b))
    C(b(i),:) = [1,0,0];
end

% [U,K,f,strain,stress,VM] = linelas3d_tetrahedron(V,T,b);
% save(saved_mats,'U', 'K', 'f', 'strain', 'stress', 'VM');


bridge = matfile(saved_mats);
U = bridge.U;
K = bridge.K;
f = bridge.f;
strain = bridge.strain;
stress = bridge.stress;
VM = bridge.VM;



% enforce dirichlet boundary
for i = 1:size(b,1)
    K(3*i-2,3*i-2) = 1.e+6;
    K(3*i-1,3*i-1) = 1.e+6;
    K(3*i,3*i) = 1.e+6;
end



% Kf = full(K);

% ErisK = readmatrix('../data/ErisK.txt','Delimiter',' ');
% ErisK = sparse(ErisK(:,1),ErisK(:,2),ErisK(:,3));
% Erisf = readmatrix('../data/Erisf.txt');

% IJVV = matdiff(K,ErisK,10);
% size(IJVV)
% fIJVV = matdiff(f,Erisf,1);
% size(fIJVV)


% ErisKf = full(ErisK);
% Kf(1:10,1:10)-ErisKf(1:10,1:10)

% f(1:10)



% options.face_vertex_color = C;
% options.face_vertex_color = U(:,2);
% options.face_vertex_color = strain(:,1).*strain(:,2).*strain(:,3);
options.face_vertex_color = VM;
plot_mesh(V,F,options);
colormap('jet');
