function [U_interp,strain,stress,VM] = compute_fields_hex(P,DV,V,Tet,C,u,r)
    U = zeros(size(u,1)/3,3);
    U(:,1) = u(1:3:end);
    U(:,2) = u(2:3:end);
    U(:,3) = u(3:3:end);

    [U_interp] = interpolate_hex_to_tet(P,U,DV,r,V);
    u_interp = zeros(size(V,2)*size(V,1),1);
    u_interp(1:3:end) = U_interp(:,1);
    u_interp(2:3:end) = U_interp(:,2);
    u_interp(3:3:end) = U_interp(:,3);

    B = zeros(6,12);
    Bs = zeros(size(Tet,1),6,12);

    for i = 1:size(Tet,1)

        Teti = Tet(i,:);

        % vertex positions
        v1 = V(Teti(1),:);
        v2 = V(Teti(2),:);
        v3 = V(Teti(3),:);
        v4 = V(Teti(4),:);

        % barycentric -> cartesian
        T = [
            1 1 1 1
            v1(1) v2(1) v3(1) v4(1)
            v1(2) v2(2) v3(2) v4(2)
            v1(3) v2(3) v3(3) v4(3)
        ];

        % volume of a tetrahedron
        tet_vol = (1/6)*det(T);

        % cartesian -> barycentric
        Tinv = inv(T)*6*tet_vol;

        % B matrix
        %       where \epsilon = B * u^e
        %       B   6x12
        %       u   12x1
        B(1,1:3:12) = Tinv(:,2);
        B(2,2:3:12) = Tinv(:,3);
        B(3,3:3:12) = Tinv(:,4);
        B(4,:) = circshift(B(1,:),1) + circshift(B(2,:),-1);
        B(5,:) = circshift(B(3,:),-1) + circshift(B(2,:),1);
        B(6,:) = circshift(B(3,:),-2) + circshift(B(1,:),2);
        B = B/(6*tet_vol);
        Bs(i,:,:)=B;
    end

    [strain,stress,vm]=per_element_fields(Tet,C,Bs,u_interp);

    N = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        N(Tet(i,:)) = N(Tet(i,:)) + 1;
    end

    VM = zeros(size(V,1),1);
    for i=1:size(Tet,1)
        VM(Tet(i,:),1) = VM(Tet(i,:),1) + vm(i)./N(Tet(i,:),1);
    end

end