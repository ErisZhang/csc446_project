clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/bb-bunny.mesh');
% [V,F,N] = readSTL('../data/origamix_rabbit.stl');

t = tsurf(F,V);

[W,BC,DV,Q] = voxelize(V,F,20,'Pad',1);
clf; hold on;
tsurf(F,V,'FaceColor',orange,falpha(1,0),fsoft);
tsurf(Q,DV,'FaceColor',blue,falpha(0.3,0.3),fsoft);
hold off;
axis equal;
axis vis3d;
rotate3d on;
set(gca,'Visible','off');


% set(gcf,'Color',0.94*[1 1 1]);
% t.EdgeColor = 'none';
% rotate3d on;
% 
% blue = [13 97 232]/255;
% orange = [255 148 0]/255;
% teal = [144 216 196]/255;
% white = [255 255 255]/255;
% bg_color = white;
% fg_color = blue;
% 
% axis equal;
% set(t,fphong,'FaceVertexCData',repmat(fg_color,size(V,1),1));
% set(t,fsoft);
% 
% l = light('Position',[0.2 -0.2 1]);
% 
% set(gca,'Visible','off');
% 
% set(gcf,'Color',bg_color);
% 
% s = add_shadow(t,l,'Color',bg_color*0.8,'BackgroundColor',bg_color,'Fade','infinite');
% 
% apply_ambient_occlusion(t,'AddLights',false,'SoftLighting',false);

