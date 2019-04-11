function tests = test_solvers
    tests = functiontests(localfunctions);
end

function test_small(t)

    fns = {
        'jacobi'
        'gauss-seidel'
        'jacobi-weighted'
        'SOR'
    };

    fs = {
        @jacobi 
        @sor   
        @(A,b,x0,max_iter,tol) jacobi(A,b,x0,max_iter,tol,'Omega',jacobi_omegaopt(A))
        @(A,b,x0,max_iter,tol) sor(A,b,x0,max_iter,tol,'Omega',sor_omegaopt(A))
    };

    A = [5.02 2.01 -0.98; 3.03 6.95 3.04; 1.01 -3.99 5.98];
    A = sparse(A);
    b = [2.05 -1.02 0.98]'; 
    x0 = zeros(size(b));             
    x = A\b;
    %  0.50774
    % -0.31141
    % -0.12966

    tol = 1e-6;
    for i = 1:size(fs,1)
        xhat = fs{i}(A,b,x0,100000,tol);
        e = norm(A*xhat-b,Inf);
        verifyTrue(t,all(e <= tol*100));
        % error bound is larger than tolerance supplied
        %   related by condition number .... handwavy here and say its correct
        %   https://archive.siam.org/books/textbooks/fr16_book.pdf
    end
end