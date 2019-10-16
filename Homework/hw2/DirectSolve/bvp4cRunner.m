function [sol] =  bvp4cRunner(params)
%bvp4cRunner Builds and runs bvp4c function using our params struct

% get the parameters
x0 = params.x0; 
xf = params.xf; 
h = params.h; 
alpha = params.alpha; 
beta = params.beta; 

% define grid 
N = (xf-x0)/h; 
xs = linspace(x0, xf, N);

% define functions needed by bvp4c
my_init = @(x) [params.infunc(x), params.infuncDer(x)]; % coordinating with my newton code 
init = bvpinit(xs, my_init);

my_bc = @(ya,yb)[ya(1)-alpha; yb(1)-beta];
my_rhs = @(~,y)[y(2); params.ffunc(y(1))];  % for our boundary value problem 

sol = bvp5c(my_rhs, my_bc, init);


end

