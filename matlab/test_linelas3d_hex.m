clear all;
format shortG;
warning('off','all');
set_project_paths();

[V,T,F] = readMESH('../data/bb-bunny.mesh');

[W,BC,DV,Q,r] = voxelize(V,F,30);

load = [0; -9.8; 0];

[VM_interp,VM,P,C,data] = linelas3d_hexahedron(W,load,r,DV,V,T);

% ori = DV(1,:);
% all_V = zeros(size(VM,1),3);
% all_F = [];

% for i = 1:size(P,1)
%     for j = 1:size(P,2)
%         for k = 1:size(P,3)
%             idx = P(i,j,k);
%             if idx ~= 0
%                 v_pos = [j-1,i-1,k-1].*r+ori;
%                 all_V(idx,:) = v_pos;
%             end
%         end
%     end
% end

% for i = 1:size(Q,1)
%     p_idx = [];
%     for j = 1:4
%         idx = floor((DV(Q(i,j),:)-ori)./r)+1;
%         p =  P(idx(2),idx(1),idx(3));
%         p_idx = [p_idx p];       
%     end
%     if any(p_idx==0)
%         continue;
%     end
%     all_F = [all_F
%              p_idx(1:3)
%              [p_idx(3) p_idx(4) p_idx(1)]];
% end

% for i = 1:size(VM,1)
%     if VM(i) > 0.25*max(VM)
%         VM(i) = 0.25*max(VM);
%     end
% end

% for i = 1:size(VM_interp,1)
%     if VM_interp(i) > 0.25*max(VM_interp)
%         VM_interp(i) = 0.25*max(VM_interp);
%     end
% end


options.face_vertex_color = log(VM_interp);
plot_mesh(V,F,options);
colormap('jet');
