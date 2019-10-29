% function to minimize
function [fval] = optfunction(eps,xp,K)
y0 = [1 sqrt(K*(-4)^2-eps)*1]; % Assign parameter values
[t,y] = ode45('shoot2',xp,y0,[],K,eps); % solve ODE
fval = abs(y(end,2)+sqrt(K*(4)^2-eps)*y(end,1)); % set parameter to minimize measure of convergence
end