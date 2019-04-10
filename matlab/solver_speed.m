clear all;
format shortG;
set_project_paths();

fns = {
    'jacobi'
    'jacobi-weighted'
    'gauss-seidel'
    '0.5-SOR'
    '1.5-SOR'
};

fs = {
    @jacobi 
    @(A,b,x0,max_iter,tol) jacobi(A,b,x0,max_iter,tol,'Omega',0.5)
    @sor   
    @(A,b,x0,max_iter,tol) sor(A,b,x0,max_iter,tol,'Omega',0.5)
    @(A,b,x0,max_iter,tol) sor(A,b,x0,max_iter,tol,'Omega',1.5)
};


% check speed
n = 200;
max_iters = 30;
x0 = zeros(n,1);
reps = 3;
tol = 1e-10;
ts = zeros(size(fs,1),reps);
for j = 1:reps
    A = sprand(n,n,n/30)+n*speye(n);
    b = rand(n,1);
    for i =1:size(fs,1)
        tic;
        fs{i}(A,b,x0,max_iters,tol);
        ts(i,j) = toc;
    end
end

sum(ts,2)./max_iters


% time(s)/iteration ( n = 200 )
% 0.0075471
% 0.011481
% 0.0030789
% 0.0087874
% 0.0087376