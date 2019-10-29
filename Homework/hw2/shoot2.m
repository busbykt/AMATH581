function rhs=shoot2(xp,y,dummy,K,eps)
rhs=[y(2);(xp^2-eps)*y(1)];