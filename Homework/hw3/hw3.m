% Kelton Busby
% AMATH 581 HW3

clear variables;

% Problem 1

xp=[-1 1]; % Define the span of the x domain
n=8; % Define the mesh size
dx=(xp(2)-xp(1))/n; % Define delta x

% A1 1st derivative center-diff periodic
e=ones(n,1);
A1=1/(2*dx)*spdiags([-e zeros(n,1) e], -1:1, n, n);

% A2 1st derivatve center-diff periodic
A2=A1;
A2(1,end)=1/(2*dx)*-1;
A2(end,1)=1/(2*dx)*1;

% A3 2nd derivative center-diff
e=ones(n,1);
A3=(1/dx^2)*spdiags([e -2*e e], -1:1, n, n);

% A4 2nd derivative center-diff periodic
A4=A3;
A4(1,n)=1/(dx^2)*1;
A4(n,1)=1/(dx^2)*1;

yp=[-1 1]; % Define the span of the y domain
m=n; % Define the mesh size in the y domain

% A5 Dx 2d center-diff 
A5 = kron(speye(n,m), A1);

% A6 Dx 2d center-diff periodic
A6 = kron(speye(n,m), A2);

% A7 Dy 2d center-diff
A7 = kron(A1, speye(n,m));

% A8 Dy 2d center diff periodic 
A8 = kron(A2, speye(n,m));

% A9 Dxx 2d center diff
A9 = kron(speye(n,m),A3);

% A10 Dxx 2d center diff periodic
A10 = kron(speye(n,m),A4);

% A11 Dyy 2d center-diff
A11 = kron(A3,speye(n,m));

%A12 Dyy 2d center-diff periodic;
A12 = kron(A4,speye(n,m));

% Laplacian no periodic
A13 = kron(speye(n,m),A3) + kron(A3,speye(n,m));

% Laplacian periodic
A14 = kron(speye(n,m),A4) + kron(A4,speye(n,m));

% -- Problem 2 --

n=10; % dimension of grid

% create asymmetric tridiagonal matrix
A = spdiags([ones(n,1)*-1.16, ones(n,1), ones(n,1)*.16],-1:1,n,n);

M = diag(A); % set predconditioner to diagonal of A
b = ones(n,1); % set b to a column of ones length n
k=0; % counter
x = zeros(n,1);


while norm(b-A*x)> 1e-8
    
    % update xk
    x = x-M.^-1*(b-A*x);
    k=k+1; % add 1 to the number of loops
end












