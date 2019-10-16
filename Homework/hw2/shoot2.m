function rhs=shoot2(xspan,y,dummy,K,eps)
rhs=[ y(2);-(K*4^2-eps)*y(1)];