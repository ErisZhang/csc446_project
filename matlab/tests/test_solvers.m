function tests = test_solvers
    tests = functiontests(localfunctions);
end

function test_small(t)
    A = [5.02 2.01 -0.98; 3.03 6.95 3.04; 1.01 -3.99 5.98];
    b = [2.05 -1.02 0.98]';              
    x = A\b;
    eps = 1e-2;
    for i = 1:size(fs,1)
        xhat = fs{i}(A,b,zeros(size(b)),1000,eps);
        verifyTrue(t,all(abs(xhat-x) <= eps));
    end
end

fs = {
    @jacobim
    @jacobie
    @gaussseidele
};

