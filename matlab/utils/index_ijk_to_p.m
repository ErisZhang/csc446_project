% a little data structure to map (i,j,k) to p
function [P,dof,B] = index_ijk_to_p(W)
    I = size(W,1);
    J = size(W,2);
    K = size(W,3);
    P = zeros(I+1,J+1,K+1);
    dof = 0;
    B = [];
    % for each grid point
    for i = 1:(I+1)
        for j = 1:(J+1)
            for k = 1:(K+1)
                % loop over its 8 neighboring cells
                is_mesh_point = 0;
                for a = -1:0
                    for b = -1:0
                        for c = -1:0
                            if (i+a)~=0 && (i+a)~=(I+1) && (j+b)~=0 && (j+b)~=(J+1) && (k+c)~=0 && (k+c)~=(K+1)
                                if W(i+a,j+b,k+c) == 1
                                    is_mesh_point = 1;
                                end
                            end
                        end

                    end
                end
                if is_mesh_point == 1
                    dof = dof + 1;
                    P(i,j,k) = dof;
                    if i == 1
                        B = [B dof]; % boundary points
                    end
                end
            end
        end
    end
end