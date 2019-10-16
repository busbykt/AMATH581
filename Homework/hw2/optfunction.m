% function to minimize

function [f,g] = optfunction(y0, eps, K)

[t,y] = ode45('shoot2',xp,y0,[],K, eps); % solve ODEs
abs(y(end,1)-0) < tol % check if solution is less than tolerance

f=y(1);
g=y(2);