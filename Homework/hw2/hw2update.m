%%

% Problem 2 - Direct Method

L=4; % set the limits of the domain
xspan=transpose(-L:.1:L); % set the span of x

% main program to solve u''(x) = f(x)
% create A

n=length(xspan)-1;
nxspan=xspan(2:n);
h = .1;

e = ones(n-1,1);
A = (1/h^2)*spdiags([-e 2+(h*nxspan).^2 -e], -1:1, n-1, n-1);

A(1,1)= 2/(3*h^2); A(1,2)= -2/(3*h^2); A(1,2)= 0;
A(end,end)=2/(3*h^2); A(end,end-1)=-2/(3*h^2); A(end,end-2)=0;


[V, D] = eig(full(A));
[d,ind] = sort(diag(D));

eigvals = d(1:5);

norm=trapz(nxspan,d(:,1).*d(:,1)); % calculate the normalization
Ds = D(ind,ind);
Vs = V(:,ind);

% save A7.dat a -ascii






