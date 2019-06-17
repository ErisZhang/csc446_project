clear all;
warning('off','all');
set_project_paths();

xs = linspace(-1,1,160);
ys = linspace(-0.1,0.1,20);
[X,Y] = meshgrid(xs,ys);
Z = zeros(size(X));
[V,F] = surf_to_mesh(X,Y,Z);
V = V(:,1:2);       % -> 2d mesh

ratio = 0.001;
load = [0; -9.8];
b = boundary_vertices(V,1,ratio);
b = reshape(2*repmat(b,2,1) - [1 0]',[],1); % vertex->node-wise
assert(size(b,1)==40);

[U,strain,stress,vm,VM] = linelas2d_triangle(V,F,b,load);

options.face_vertex_color = VM;
plot_mesh(V,F,options);
colormap('jet');

% tsurf(F,V);
% axis equal;



% % Fit to half the unit square
% V = V/(2*max(max(V)-min(V))); 
% % Initialize as stretched object
% U = 1.5*V-V;
% Ud = zeros(size(V));
% 
% % https://github.com/alecjacobson/gptoolbox/blob/master/mesh/tsurf.m
% %   wraps  https://www.mathworks.com/help/matlab/ref/trisurf.html
% t = tsurf(F,V+U);
% axis equal;
% axis manual;
% while true
%     %[U,Ud] = linelas2d_gptoolbox(V,F,[],[],'U0',U,'Ud0',Ud);
%     [U,Ud] = linelas2d_gptoolbox_prev(V,F,[],'U0',U,'Ud0',Ud);
% 
%     t.Vertices = V+U;
%     drawnow;
% end

