function [G, J] = buildSystemNL(theta, params)
%buildSystemNL builds nonlinear equation G and its jacobian J

% load stuff from params
ffunc = params.ffunc; % function
ffuncDer = params.ffuncDer; 
h = params.h;  % step size
x0 = params.x0; % initial domain point 
xf = params.xf; % final domain point; 
% initial conditions
alpha = params.alpha;
beta = params.beta; 


N = (xf-x0)/h; % number of steps 
assert(N-round(N) == 0); % make sure number of steps is even. 



%sparse matrix for finite differences
e = ones(N,1);
A = (1/h^2)*spdiags([e -2*e e], -1:1, N, N); 
G = A*theta - ffunc(theta);
G(1) = G(1) + alpha/h^2;
G(end) = G(end) + beta/h^2;


J = (1/h^2)*spdiags([e (-2*e - h^2*ffuncDer(theta)) e], -1:1, N, N); 






end

