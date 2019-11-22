% Kelton Busby
% AMATH 581 HW4

clear variables;

% set value of kappa
kappa=.1;

% set space and time discretization
dt = .05;
dx = .1;

tspan = 0:dt:1;
xspan = 0:dx:1;
n = length(tspan);
m = length(xspan);


% initial distribution
init_cond = @(x) sin(2*pi*x);

% define initial U
U_fe = transpose(init_cond(xspan));

r = kappa*dt/(dx^2);

e=ones(length(xspan),1);
A = spdiags([e -2*e e], -1:1, m, m);
A(1,10)=1;
A(end,2)=1;

U = zeros(length(xspan),length(tspan));

for i=1:length(tspan)
    U(:,i)=U_fe;
    U_fe = U_fe + r*A*U_fe;
    plot(U)
    hold on
    drawnow;
    pause(.1);
end

% save the output of final time t=1 A1
A1 = U(:,end);
save A1.dat A1 -ascii



%%

% backward Euler Method
% set space and time discretization
dt = .05;
dx = .1;
U_be = transpose(init_cond(xspan)); % define initial U

U = zeros(length(xspan),length(tspan));

for i=1:length(tspan)
    
    M = diag(diag(A)); % set predconditioner to diagonal of A
    b = ones(m,1); % set b to a column of ones length n
    x = U_be; % set initial guess
    k=1; % counter
    errVec = zeros(1,1000); % initialize array to hold errors
    resid = x - r*A*x - U_be; % residuals
    err = norm(resid); % compute the error from the residuals
    errVec(k)=err;
    z = M\resid; % solve for z
    
    while err > 1e-8
        x = x + z; % update x
        resid = x - r*A*x - U_be; % update residuals
        err = norm(resid); % update error
        z = M\resid; % update z
        k=k+1; % add 1 to the number of loops
        errVec(k)=err;
        if k > 1000
            break
        end
    end
    U(:,i)=U_be;
    U_be = x; % update U_be
    plot(U)
    hold on
    drawnow;
    pause(.1);
end

%Save A2
A2 =U(:,end);
save A2.dat A2 -ascii


%%

% backward Euler Method small step
% set space and time discretization
dt = .01;
dx = .01;
tspan = 0:dt:1;
xspan = 0:dx:1;
n = length(tspan);
m = length(xspan);

r = kappa*dt/(dx^2);

e=ones(length(xspan),1);
A = spdiags([e -2*e e], -1:1, m, m);
A(1,m-1)=1;
A(end,2)=1;

U_be = transpose(init_cond(xspan)); % define initial U
U = zeros(length(xspan),length(tspan));

for i=1:length(tspan)
    
    M = diag(diag(A)); % set predconditioner to diagonal of A
    b = ones(m,1); % set b to a column of ones length n
    x = U_be; % set initial guess
    k=1; % counter
    errVec = zeros(1,1000); % initialize array to hold errors
    resid =x - r*A*x - U_be; % residuals
    err = norm(resid); % compute the error from the residuals
    errVec(k)=err;
    z = M\resid; % solve for z
    
    while err > 1e-8
        x = x + z; % update x
        resid = x - r*A*x - U_be; % update residuals
        err = norm(resid); % update error
        z = M\resid; % update z
        k=k+1; % add 1 to the number of loops
        errVec(k)=err;
        if k > 1000
            break
        end
    end
    U(:,i)=U_be;
    U_be = x; % update U_be
    plot(U)
    hold on
    drawnow;
    pause(.01);
end

%%

% Crank-Nicolson Method
% set space and time discretization
dt = .05;
dx = .1;
tspan = 0:dt:1;
xspan = 0:dx:1;
n = length(tspan);
m = length(xspan);

r = kappa*dt/(dx^2);

e=ones(length(xspan),1);
A = spdiags([e -2*e e], -1:1, m, m);
A(1,m-1)=1;
A(end,2)=1;
U_cn = transpose(init_cond(xspan)); % define initial U

U = zeros(length(xspan),length(tspan));

for i=1:length(tspan)
    
    M = diag(diag(A)); % set predconditioner to diagonal of A
    b = ones(m,1); % set b to a column of ones length n
    x = U_cn; % set initial guess
    k=1; % counter
    errVec = zeros(1,1000); % initialize array to hold errors
    resid = x - r*A*(x+U_cn) - U_cn; % residuals
    err = norm(resid); % compute the error from the residuals
    errVec(k)=err;
    z = M\resid; % solve for z
    
    while err > 1e-8
        x = x + z; % update x
        resid = x - r*A*(x+U_cn) - U_cn; % update residuals
        err = norm(resid); % update error
        z = M\resid; % update z
        k=k+1; % add 1 to the number of loops
        errVec(k)=err;
        if k > 1000
            break
        end
    end
    U(:,i)=U_cn;
    U_cn = x; % update U_be
    plot(U)
    hold on
    drawnow;
    pause(.1);
end

%Save A4
A4 = U(:,end);
save A4.dat A4 -ascii

%%
clear variables;
% Crank-Nicolson Method
% set space and time discretization
% initial distribution
init_cond = @(x) sin(2*pi*x);
kappa=.1;
dt = .01;
dx = .01;
tspan = 0:dt:1;
xspan = 0:dx:1;
n = length(tspan);
m = length(xspan);

r = kappa*dt/(dx^2);
U_cn = transpose(init_cond(xspan)); % define initial U

e=ones(length(xspan),1);
A = spdiags([e -2*e e], -1:1, m, m);
A(1,end-1)=1;
A(end,2)=1;

U = zeros(length(xspan),length(tspan));

for i=1:length(tspan)
    
    M = diag(diag(A)); % set predconditioner to diagonal of A
    b = ones(m,1); % set b to a column of ones length n
    x = U_cn; % set initial guess
    k=1; % counter
    errVec = zeros(1,1000); % initialize array to hold errors
    resid = x - r*A*(x+U_cn) - U_cn; % residuals
    err = norm(resid); % compute the error from the residuals
    errVec(k)=err;
    z = M\resid; % solve for z
    
    while err > 1e-2
        x = x + z; % update x
        resid = x - r*A*(x+U_cn) - U_cn; % update residuals
        err = norm(resid); % update error
        z = M\resid; % update z
        k=k+1; % add 1 to the number of loops
        errVec(k)=err;
        if k > 1000
            break
        end
    end
    U(:,i)=U_cn;
    U_cn = x; % update U_be
    plot(U)
    hold on
    drawnow;
    pause(.1);
end

