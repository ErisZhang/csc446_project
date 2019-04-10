clear all;
format shortG;
set_project_paths();


n = 50;
A = rand(n,n)+n*eye(n);
b = rand(n,1);
x = A\b;

saved_mats = '../data/archbridge_withbc.mat';
bridge = matfile(saved_mats);
U = bridge.U;
K = bridge.K;
f = bridge.f;
u = K\f;


x0 = zeros(size(f));
max_iters = 2000;
tol = 1e-14;

uhat = sor(K,f,x0,max_iters,tol,1.8);
IJVV = matdiff(u,uhat,1e-6)
size(IJVV)



