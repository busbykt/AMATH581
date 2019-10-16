function [A, F] = buildSystem(params)
%buildSystem builds sparse matrix A and right hand side F

% load stuff from params
ffunc = params.ffunc; % function for right hand side 
h = params.h;  % step size
x0 = params.x0; % initial domain point 
xf = params.xf; % final domain point; 
% initial conditions
alpha = params.alpha;
beta = params.beta; 


N = (xf-x0)/h; % number of steps 
assert(N-round(N) == 0); % make sure number of steps is an integer. 



xs = linspace(x0, xf, N); 
F = ffunc(xs)'; 
F(1) = F(1) - alpha/h^2;
F(end) = F(end) - beta/h^2;

%fancy way to build sparse matrix
e = ones(N,1);
A = (1/h^2)*spdiags([e -2*e e], -1:1, N, N); 








end

