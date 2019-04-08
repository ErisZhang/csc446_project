% find index and corresponding elements that are
%       different when comparing `A` and `B`
%       up to tolerence `tol`
function IJVV = matdiff(A,B,tol)
    assert(all(size(A) == size(B)), 'A,B should be same size');

    IJVV = zeros(0,4);
    ndiff = 1;

    if issparse(A) && issparse(B)
        [AI,AJ,AV] = find(A);
        AIJV = [AI,AJ,AV];
        [BI,BJ,BV] = find(B);
        BIJV = [BI,BJ,BV];
        
        for k = 1:size(AIJV,1)
            if mod(k,100000) == 0
                [k size(AIJV,1)]
            end
            i = AIJV(k,1);
            j = AIJV(k,2);
            v = AIJV(k,3);
            Bij = full(B(i,j));
            if abs(Bij-v)>tol
                IJVV(ndiff,:) = [i,j,v,Bij];
                ndiff = ndiff + 1;
            end
        end

        for k = 1:size(BIJV,1)
            if mod(k,1000) == 0
                [k size(BIJV,1)]
            end
            i = BIJV(k,1);
            j = BIJV(k,2);
            v = BIJV(k,3);
            Aij = full(A(i,j));
            if abs(Aij-v)>tol
                IJVV(ndiff,:) = [i,j,Aij,v];
                ndiff = ndiff + 1;
            end
        end
    else
        for i = 1:size(A,1)
            for j = 1:size(A,2)
                if abs(A(i,j)-B(i,j))>tol
                    IJVV(ndiff,:) = [i j A(i,j) B(i,j)];
                    ndiff = ndiff + 1;
                end
            end
        end
    end
end