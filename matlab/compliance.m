% import mesh

[V,F] = bwmesh('Documents/MATLAB/hidden-supports-2d/example-scene-0.1x.png');
%[V,F] = bwmesh('Documents/MATLAB/hidden-supports-2d/,...
%   example-scene-supported-0.1x.png');

% create grid
num_voxels = 70;
[G,side,r] = voxel_grid(V,num_voxels);

% generate views
num_views = 10;
max_y = 90;
min_y = 45;
mu = (min_y+max_y)/2;
sigma = abs((max_y-min_y)/2);
views = normrnd(mu,sigma,[num_views, 1]);

% intersect rays with mesh
E = edges(F); % get edges

% place to hold visibility values
S = zeros(size(G,1),1);

for view = 1:num_views
    o = [140 views(view)];
    view
    o
    %hold on
    %axis([0 148 0 140])
    %plot(o(1),o(2),'r*');
    %drawnow;

    for i = 1:size(G,1)
        dir = G(i,:)-o;
        dist = sqrt(sum(dir.^2,2));
        dir = normalizerow(dir);
        [intersect, t, lambda] = ray_polygon_intersect(o,dir,V,E);
        if any(t(intersect)<dist)
            S(i) = S(i) + 1;
        end
    end
end
S = S / num_views;


%%
V_list = cell(10,1);
F_list = cell(10,1);
stress_list = cell(10,1);
for i=1:10
    i
    S_i = S;
    S_i(S<(i/10) & S>0) = 0;
    
    % marching squares
    s = reshape(S_i,num_voxels,num_voxels,1);
    s = flip(s);
    [V_s,F_s] = bwmesh(s);

    % scale   
    scale = max(max(V)) / max(max(V_s));
    V_s = scale .* V_s;

    x = V_s(:,1);
    y = V_s(:,2);

    % get boundaries
    % this will be a logical matrix into V_s
    b_bottom = (y >= 0.5 & y <= 6.5);
    b_left = (x >= 0.5 & x <= 6.5);
    b = b_bottom | b_left;

    % matrix of body forces (gravity)
    % -9.8 m/s�
    bf = [zeros(size(V_s,1),1),(-9.8)*ones(size(V_s,1),1)];

    % initial run of linear elasticity
    % Young's modulus 2.18e+9 Pa
    [U, Ud, data] = linear_elasticity_gptoolbox_2d(V_s,F_s,b, ...
        'BodyForces',bf,'Young',2.18e+9,'Nu',0.35);
    
    V_list{i} = V_s;
    F_list{i} = F_s;
    % for ABS yield strength is 2.6e+7 Pa to 4.48e+7 Pa
    % if the stress under 
    stress_list{i} = data.VM;

end

%%
% visualize
%tsh = tsurf(F_s,V_s,'CDATA',normrow(U));
%tsh = tsurf(F,V);

% get max stress
max_stress = cellfun(@max,stress_list);

subplot(1,2,1);
tsurf(F_list{1},V_list{1},'CDATA',stress_list{1});
caxis([0 max(max_stress)]);
axis equal;

subplot(1,2,2);
tsurf(F_list{10},V_list{10},'CDATA',stress_list{10});
max_stress
caxis([0 max(max_stress)]);
axis equal;

