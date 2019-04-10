clear all;
format shortG;
set_project_paths();

fs = {
    @jacobim
    @jacobie
    @(A,b,x0,max_iters,tol) jacobis(A,b,x0,max_iters,tol,2/3)
    @gaussseidele
    @gaussseidels
};


% check speed
n = 200;
max_iters = 30;
x0 = zeros(n,1);
reps = 3;
ts = zeros(size(fs,1),reps);
for j = 1:reps
    A = rand(n,n)+n*eye(n);
    b = rand(n,1);
    for i =1:size(fs,1)
        tic;
        fs{i}(A,b,x0,max_iters,tol);
        ts(i,j) = toc;
    end
end

sum(ts,2)./max_iters
%
% time(s)/iteration ( n = 200 )
% 0.000146            jacobim
% 0.00025663          jacobie
% 0.00023259          gaussseidele