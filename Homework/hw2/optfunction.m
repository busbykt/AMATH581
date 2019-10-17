% function to minimize

function f = optfunction(y0, eps)

[t,y] = ode45('shoot2',xp,y0,[],1, eps); % solve ODEs
f = abs(y(end,1)-0); % set parameter to minimize