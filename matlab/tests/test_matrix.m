function tests = test_isspd
    tests = functiontests(localfunctions);
end

function test_eye(t)
    A = eye(2);
    verifyEqual(t,isspd(A),true)
end

function test_real(t)
    A = [2 -1 0; -1 2 -1; 0 -1 2];
    verifyEqual(t, isspd(A), true);
end

function test_diagonally_dominant(t)
    v = zeros(10,1);
    B = [v-1, v+4,v-1];
    d = [-1,0,1];
    S = spdiags(B,d,10,10);
    verifyEqual(t, isspd(full(S)), true);
end

function test_not_diagonally_dominant(t)
    v = zeros(10,1);
    B = [v+2,v+2,v+2];
    d = [-1,0,1];
    S = spdiags(B,d,10,10);
    verifyEqual(t, isspd(full(S)), false);
end
