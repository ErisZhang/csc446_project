function tests = test_matrix
    tests = functiontests(localfunctions);
end

function test_isspd_eye(t)
    A = eye(2);
    verifyEqual(t,isspd(A),true)
end

function test_isspd_real(t)
    A = [2 -1 0; -1 2 -1; 0 -1 2];
    verifyEqual(t, true,isspd(A));
end

function test_isspd_diagonally_dominant(t)
    v = zeros(10,1);
    B = [v-1, v+4,v-1];
    d = [-1,0,1];
    S = spdiags(B,d,10,10);
    verifyEqual(t, true, isspd(full(S)));
end

function test_isspd_not_diagonally_dominant(t)
    v = zeros(10,1);
    B = [v+2,v+2,v+2];
    d = [-1,0,1];
    S = spdiags(B,d,10,10);
    verifyEqual(t, false, isspd(full(S)));
end

function test_isdiagonallydonimantrow_true(t)
    A = [
        -2 1 -0.5
        1 2 -0.9
        0 0 1
    ];
    [b,l]=isdiagdominbyrow(A);
    verifyEqual(t,true,b);
end

function test_isdiagonallydonimantrow_false(t)
    A = [
        -1.5 1 -0.5
        1 2 -0.9
        0 0 1
    ];
    [b,l]=isdiagdominbyrow(A);
    verifyEqual(t,false,b);
end

function test_nnz_indices_byrow(t)
    n = 10;
    S = sprand(n,n,n/30);
    indices = nnz_indices_byrow(S);
    for i = size(indices,1)
        for j = size(indices{i},2)
            verifyTrue(t, S(i,indices{i}(j)) ~= 0);
        end
    end
end