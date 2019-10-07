function dydt = vanderpol(t,y,eps)
%VDP1  Evaluate the van der Pol ODEs for eps = 1

dydt = [y(2); (1-y(1)^2)*y(2)-y(1)];
