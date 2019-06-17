function [U,strain,stress,vm,VM] = compute_fields_2d(V,F,D,Bs,u)
    U = zeros(size(V));
    U(:,1) = u(1:2:end);
    U(:,2) = u(2:2:end);

%     [strain,stress,vm]=per_element_fields(F,D,Bs,u);
    
    ij2p = zeros(6,1); % for indexing from `V` to `K/f/u`
    strain = zeros(size(F,1),3);
    stress = zeros(size(F,1),3);
    vm     = zeros(size(F,1),1);
    
    for i = 1:size(F,1)
        B = squeeze(Bs(i,:,:)); % with dimensions of length 1 removed
        ele_i = F(i,:);
        ij2p(1:2:end) = 2*ele_i-1;
        ij2p(2:2:end) = 2*ele_i-0;
        strain(i,:) = B*u(ij2p);
        stress(i,:) = D*strain(i,:)';
        
        vm(i) = sqrt(((stress(i,1)-stress(i,2))^2+(stress(i,2)-stress(i,3))^2+(stress(i,3)-stress(i,1))^2)/2);
    end
    
   
    N = zeros(size(V,1),1);
    
    for i=1:size(F,1)
        N(F(i,:)) = N(F(i,:)) + 1;
    end
    
    VM = zeros(size(V,1),1);
    
    for i=1:size(F,1)
        VM(F(i,:),1) = VM(F(i,:),1) + vm(i)./N(F(i,:),1);
    end
end
