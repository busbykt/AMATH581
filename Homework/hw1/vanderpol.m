function dydt = vanderpol(t,y, eps)
%Evaluate the van der Pol ODEs for provided eps

dydt = [y(2); eps*(1-y(1)^2)*y(2)-y(1)];
