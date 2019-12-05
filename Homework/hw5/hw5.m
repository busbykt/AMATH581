% Kelton Busby
% AMATH 581 HW5

clear variables;

% Problem 1

% define the span of x
xspan = -1:.05:1;
% define some initial guess for the solution
init = bvpinit(xspan,[0 0]);
% solve the boundary value problem
sol = bvp4c(@bvp_rhs, @bvp_bc, init);

x = linspace(-1,1,100); BS=deval(sol,x);
plot(x,BS(1,:))

%save the solution at u(0)
A1 = deval(sol,0);
A1 = A1(1);
save A1.dat A1 -ascii

%%

% Problem 2

N=16;
[D,x] = cheb(N); D2 = D*D; D2 = D2(2:N,2:N);
u = zeros(N-1,1);
I=eye(length(D2));
change = 1; it = 0;
while change > 1e-15 % fixed point iteration
    J = D2-exp(u); % what is the jacobian??
    unew = u - J^-1*exp(u);
    change = norm(unew-u,inf);
    u = unew; it=it+1;
end
u = [0;u;0];
clf, subplot('position',[.1 .4 .8 .5])
plot(x,u,'.','markersize',16)
xx = -1:.01:1;
uu = polyval(polyfit(x,u,N),xx);
line(xx,uu), grid on
title(sprintf('no steps = %d   u(0)=%18.14f',it,u(N/2+1)))

% Problem 3

% Need to convert a PDE to an ODE?