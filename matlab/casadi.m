addpath('/Users/zhangjiayi/Documents/github/casadi');
import casadi.*;

x = SX.sym('x');
obj = x^2-6*x+13; % symbolic expression

g = []; % optimization constraints - empty
P = []; % optimization problem parameters - empty

OPT_variables = x;
nlp_prob = struct('f',obj,'x',OPT_variables,'g',g,'p',P);

opts = struct;
opts.ipot.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

% Ipopt: interior Point Optimizer
solver = nlpsol('solver','ipopt',nlp_prob,opts);



% a = SX.sym('a',2);
% b = SX.sym('b',2);
% 
% states = [a,b];
% speed = norm_2(a)+norm_2(b);
% 
% rhs = [a;speed];
% J = jacobian(rhs,states);
% f = MXFunction('f', {states}, {rhs});
% f({[1;1;1;1]})

% V = 4/3*pi*r^3;
% A = jacobian(V,r);
% f = Function('f', {r}, {A});
% f(1);