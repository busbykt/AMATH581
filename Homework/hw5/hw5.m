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
%plot(x,BS(1,:))

%save the solution at u(0)
A1 = deval(sol,0);
A1 = A1(1);
save A1.dat A1 -ascii


% Problem 2

N=16;
[D,x] = cheb(N); D2 = D*D; D2 = D2(2:N,2:N);
u = zeros(N-1,1);
Ix=eye(length(D2));
change = 1; it = 0;
save_changes=zeros(100,1);
while change > 1e-15 % fixed point iteration
    J = D2 - diag(exp(u)); % Define the Jacobian
    f = D2*u-exp(u);
    unew = u-J\f;
    change = norm(unew-u,inf);
    u = unew; it=it+1;
    save_changes(it) = change;
end

%save the changes at iteration 1 and 2
A2 = save_changes(1);
A3 = save_changes(2);
save A2.dat A2 -ascii
save A3.dat A3 -ascii

% Problem 3

clear variables;
tspan = 0:.0001:3.55;
xspan = -1:.02:1;
N = length(xspan);
[D,x] = cheb(N); D2 = D*D; D2 = D2(2:N,2:N);
u0 = zeros(N-1,1);

[t,u] = ode23s(@(t,u)D2*u+exp(u),tspan,u0);

% get the index where x = 0 for A4
Ix = find(xspan == 0);

% get the index where t = 3.5 for A4
It = find(tspan == 3.5);

A4 = u(It,Ix);
save A4.dat A4 -ascii

% get the values before and after 5
t0 = find(u(:,Ix) > 4.99 & u(:,Ix) < 5);
t1 = find(u(:,Ix) > 5 & u(:,Ix) < 5.01);

%slope
slope = (u(t1,Ix)-u(t0,Ix))/((t1-t0)*(tspan(end)-tspan(1))/length(tspan));

dt = (5-u(t0,Ix))/slope;

tval = tspan(t0)+dt;

A5 = tval;
save A5.dat A5 -ascii

