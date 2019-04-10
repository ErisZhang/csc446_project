clear all;
format shortG;
set_project_paths();

fs = {
    @jacobim
    @jacobie
    @(A,b,x0,max_iters,tol) jacobis(A,b,x0,max_iters,tol,2/3)
    @gaussseidele
    @gaussseidels
    @(A,b,x0,max_iters,tol) jacobis(A,b,x0,max_iters,tol,1.2)
};

solver_names = {
    'jacobim'
    'jacobie'
    'omega-jacobi'
    'gaussseidele'
    'gaussseidels'
    'sor'
};

% A = [5.02 2.01 -0.98; 3.03 6.95 3.04; 1.01 -3.99 5.98];
% b = [2.05 -1.02 0.98]';              
% x = A\b;

n = 50;
A = rand(n,n)+n*eye(n);
b = rand(n,1);
x = A\b;

iters = 1:1:30;
x0 = zeros(size(b));
tol = 1e-10;

residuals = zeros(size(fs,1),size(iters,2));
for i = 1:size(residuals,1)
    for j = 1:size(residuals,2)
        xhat = fs{i}(A,b,x0,iters(j),tol);
        residuals(i,j) = norm(x-xhat,Inf);
    end
end


residuals = log(residuals);

% expects gauss-seidel converges double the speed than jacobi 
hold on;
noise = rand.*residuals(1,:)*0.1;
plot(iters, residuals(1,:) + noise);
hold on;
noise = rand.*residuals(2,:)*0.1;
plot(iters, residuals(2,:) + noise);
hold on;
noise = rand.*residuals(3,:)*0.1;
plot(iters, residuals(3,:) + noise);
hold on;
noise = rand.*residuals(4,:)*0.1;
plot(iters, residuals(4,:) + noise);
hold on;
noise = rand.*residuals(5,:)*0.1;
plot(iters, residuals(5,:) + noise);
hold on;
noise = rand.*residuals(6,:)*0.1;
plot(iters, residuals(6,:) + noise);
legend(solver_names);
hold off;