function [U_interp] = interpolate_hex_to_tet(P,U,DV,r,V)
% myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description
    U_interp = zeros(size(V,1),1);

    ori = DV(1,:);
    idx = floor((V-ori)./r)+1;
    
    cell_center = (idx-0.5).*r+ori;
    local_coor = (V-cell_center)./(r/2);

    for m = 1:size(V,1)
        a = local_coor(1);
        b = local_coor(2);
        c = local_coor(3);

        u_interp = [0];
        
        for i = 0:2:2
            for j = 0:2:2
                for k = 0:2:2
                    index = idx(m,:)+[i/2 j/2 k/2];
                    P_index = P(index(2),index(1),index(3));
                    displacement = U(P_index,:);
                    weight = 1/8*(1+(i-1)*a)*(1+(j-1)*b)*(1+(k-1)*c);
                    u_interp = u_interp+displacement*weight;
                end
            end
        end

        U_interp(m,:) = u_interp;

    end
end