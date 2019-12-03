% Kelton Busby
% AMATH 581 HW5

clear variables;

% Problem 1

% define the span of x
xspan = -1:.1:1;
% define some initial guess for the solution
init = bvpinit(xspan,[0 0]);
% solve the boundary value problem
sol = bvp4c(@bvp_rhs, @bvp_bc, init);

x = linspace(-1,1,100); BS=deval(sol,x);
plot(x,BS(1,:))


% Problem 2

