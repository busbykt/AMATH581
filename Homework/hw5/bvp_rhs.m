%bvp_rhs
function rhs=bvp_rhs(x,y)
rhs=[y(2); sin(8*x)-4*y(2)-exp(x)*y(1)];