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
A1s = full(A1);
save A1.dat A1s -ascii

% A2 1st derivatve center-diff periodic
A2=A1;
A2(1,end)=1/(2*dx)*-1;
A2(end,1)=1/(2*dx)*1;
A2s = full(A2);
save A2.dat A2s -ascii

% A3 2nd derivative center-diff
e=ones(n,1);
A3=(1/dx^2)*spdiags([e -2*e e], -1:1, n, n);
A3s = full(A3);
save A3.dat A3s -ascii

% A4 2nd derivative center-diff periodic
A4=A3;
A4(1,n)=1/(dx^2)*1;
A4(n,1)=1/(dx^2)*1;
A4s = full(A4);
save A4.dat A4s -ascii

yp=[-1 1]; % Define the span of the y domain
m=n; % Define the mesh size in the y domain

% A5 Dx 2d center-diff 
A5 = kron(speye(n,m), A1);
A5s = full(A5);
save A5.dat A5s -ascii

% A6 Dx 2d center-diff periodic
A6 = kron(speye(n,m), A2);
A6s = full(A6);
save A6.dat A6s -ascii

% A7 Dy 2d center-diff
A7 = kron(A1, speye(n,m));
A7s = full(A7);
save A7.dat A7s -ascii

% A8 Dy 2d center diff periodic 
A8 = kron(A2, speye(n,m));
A8s = full(A8);
save A8.dat A8s -ascii

% A9 Dxx 2d center diff
A9 = kron(speye(n,m),A3);
A9s = full(A9);
save A9.dat A9s -ascii

% A10 Dxx 2d center diff periodic
A10 = kron(speye(n,m),A4);
A10s = full(A10);
save A10.dat A10s -ascii

% A11 Dyy 2d center-diff
A11 = kron(A3,speye(n,m));
A11s = full(A11);
save A11.dat A11s -ascii

%A12 Dyy 2d center-diff periodic;
A12 = kron(A4,speye(n,m));
A12s = full(A12);
save A12.dat A12s -ascii

% Laplacian no periodic
A13 = kron(speye(n,m),A3) + kron(A3,speye(n,m));
A13s = full(A13);
save A13.dat A13s -ascii

% Laplacian periodic
A14 = kron(speye(n,m),A4) + kron(A4,speye(n,m));
A14s = full(A14);
save A14.dat A14s -ascii

% -- Problem 2 --

n=10; % dimension of grid

% create asymmetric tridiagonal matrix
A = spdiags([ones(n,1)*-1.16, ones(n,1), ones(n,1)*.16],-1:1,n,n);

M = diag(diag(A)); % set predconditioner to diagonal of A
b = ones(n,1); % set b to a column of ones length n
k=1; % counter
x = zeros(n,1); % set initial guess
errVec = zeros(1,1000); % initialize array to hold errors
r = b-A*x; % residuals
err = norm(r); % compute the error from the residuals
errVec(k)=err;
z = M\r; % solve for z

while err > 1e-8
    x = x + z; % update x
    r = b - A*x; % update residuals
    err = norm(r); % update error
    z = M\r; % update z
    k=k+1; % add 1 to the number of loops
    errVec(k)=err;
    if k > 1000
        break
    end
end

errVec = errVec(1:k);

% output answers to file
save A15.dat errVec -ascii

n=30; % dimension of grid

% create asymmetric tridiagonal matrix
A = spdiags([ones(n,1)*-1.16, ones(n,1), ones(n,1)*.16],-1:1,n,n);

M = diag(diag(A)); % set predconditioner to diagonal of A
b = ones(n,1); % set b to a column of ones length n
k=1; % counter
x = zeros(n,1); % set initial guess
errVec30 = zeros(1,1000); % initialize array to hold errors
r = b-A*x; % residuals
err = norm(r); % compute the error from the residuals
errVec30(k)=err;
z = M\r; % solve for z

while err > 1e-8
    x = x + z; % update x
    r = b - A*x; % update residuals
    err = norm(r); % update error
    z = M\r; % update z
    k=k+1; % add 1 to the number of loops
    errVec30(k)=err;
    if k > 1000
        break
    end
end

errVec30 = errVec30(1:k);

% output answers to file
save A16.dat errVec30 -ascii

% Gauss-Seidel
n=10; % dimension of grid

% create asymmetric tridiagonal matrix
A = spdiags([ones(n,1)*-1.16, ones(n,1), ones(n,1)*.16],-1:1,n,n);

M = tril(A); % set predconditioner to diagonal of A
b = ones(n,1); % set b to a column of ones length n
k=1; % counter
x = zeros(n,1); % set initial guess
gs10errVec = zeros(1,1000); % initialize array to hold errors
r = b-A*x; % residuals
err = norm(r); % compute the error from the residuals
gs10errVec(k)=err;
z = M\r; % solve for z

while err > 1e-8
    x = x + z; % update x
    r = b - A*x; % update residuals
    err = norm(r); % update error
    z = M\r; % update z
    k=k+1; % add 1 to the number of loops
    gs10errVec(k)=err;
    if k > 1000
        break
    end
end

gs10errVec = gs10errVec(1:k);

% output answers to file
save A17.dat gs10errVec -ascii

n=30; % dimension of grid

% create asymmetric tridiagonal matrix
A = spdiags([ones(n,1)*-1.16, ones(n,1), ones(n,1)*.16],-1:1,n,n);

M = tril(A); % set predconditioner to diagonal of A
b = ones(n,1); % set b to a column of ones length n
k=1; % counter
x = zeros(n,1); % set initial guess
gs30errVec = zeros(1,1000); % initialize array to hold errors
r = b-A*x; % residuals
err = norm(r); % compute the error from the residuals
gs30errVec(k)=err;
z = M\r; % solve for z

while err > 1e-8
    x = x + z; % update x
    r = b - A*x; % update residuals
    err = norm(r); % update error
    z = M\r; % update z
    k=k+1; % add 1 to the number of loops
    gs30errVec(k)=err;
    if k > 1000
        break
    end
end

gs30errVec = gs30errVec(1:k);

% output answers to file
save A18.dat gs30errVec -ascii


% plot(errVec)
% hold on
% plot(gs10errVec)
% 
% figure
% plot(errVec30)
% hold on
% plot(gs30errVec)

