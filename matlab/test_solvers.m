clear all;
format shortG;
set_project_paths();

fs = {
    @jacobim
    @jacobie
    @gaussseidele
};

A = [5.02 2.01 -0.98; 3.03 6.95 3.04; 1.01 -3.99 5.98];
b = [2.05 -1.02 0.98]';              
x = A\b;


iters = 1:1:30
x0 = zeros(size(b));
eps = 1e-20;

% solvers x iterations x residual
residuals = zeros(size(fs,1),size(iters,2));

for i = 1:size(residuals,1)
    for j = 1:size(residuals,2)
        xhat = fs{i}(A,b,x0,iters(j),eps);
        residuals(i,j) = norm(abs(x-xhat),Inf);
    end
end

residuals = log(residuals);

plot(iters,residuals(1,:));
hold on;
plot(iters,residuals(2,:));
hold on;
plot(iters,residuals(3,:));

legend('jacobim', 'jacobie', 'gaussseidele');

% expects gauss-seidel converges double the speed than jacobi 