% HW2
% Question 1

%Initialization
clear all;
close all;

% set a tolerance for shooting
tol = 10^(-4);

% set eigenfunction colors
col = ['r','b','g','c','m'];


K=1; % defining parameter K to keep integral equal to 1?
L=4; % set boundary
xp=-4:.1:4; % span of the computational domain

eigenvalues = zeros(5,1);
eigenvectors = zeros(81,5);

eps_start=0; % beginning value for epsilon
for modes=1:5 % begin mode loop
    eps=eps_start; % initial value of the eigenvalue epsilon
    deps=1; % default step size in epsilon - set by Sasha
    for j=1:1000 % begin convergence loop for epsilon
        y0=[1 sqrt(K*L^2-eps)*1]; % initial conditions
        [t,y] = ode45('shoot2',xp,y0,[],K,eps); % solve ODE
        err = y(end,2)+sqrt(K*L^2-eps)*y(end,1); % calculate error
        %err = y(end,1);
        if abs(err) < tol
            eigenvalues(modes,1)=eps; % save the eigenvalue
            break % leave convergence loop
        end
        if (-1)^(modes+1)*err>0 % check if epsilon should be higher or lower
            eps=eps+deps; % set epsilon lower by delta epsilon
            deps=deps/2; % set the change in epsilon to 1/2 previous change
        else
            eps=eps-deps; % set epsilon higher by delta epsilon
        end
    end % end converge loop
    eps_start=eps+2; % after finding eigenvalue pick new starting value for next mode
    norm=trapz(t,y(:,1).*y(:,1)); % calculate the normalization 
    %plot(t,y(:,1)/sqrt(norm),col(modes)); hold on % plot each mode
    eigenvectors(:,modes)=abs(y(:,1)/sqrt(norm)); % save the eigenvector
end
%legend('mode1','mode2','mode3','mode4','mode5');

% output answers to file
a = eigenvectors(:,1);
save A1.dat a -ascii
a = eigenvectors(:,2);
save A2.dat a -ascii
a = eigenvectors(:,3);
save A3.dat a -ascii
a = eigenvectors(:,4);
save A4.dat a -ascii
a = eigenvectors(:,5);
save A5.dat a -ascii
save A6.dat eigenvalues -ascii

% Problem 2

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

first_row=zeros(1,5);
last_row=zeros(1,5);

for modes=1:5
    first=(2*Vs(1,modes)-.5*Vs(2,modes))/(h*sqrt(4^2-Ds(modes,modes))+3/2);
    last=(-2*Vs(end,end)+.5*Vs(modes,2))/(h*sqrt(4^2-Ds(modes,modes))-3/2);
    
    first_row(modes)=first;
    last_row(modes)=last;
    
end

first_row = [first_row,zeros(1,n-1-modes)];  % horizontal concatenation
last_row = [last_row,zeros(1,n-1-modes)];  % horizontal concatenation

Vs = abs([first_row;Vs;last_row]);

% output answers to file
a = eigenvectors(:,1);
save A7.dat a -ascii
a = eigenvectors(:,2);
save A8.dat a -ascii
a = eigenvectors(:,3);
save A9.dat a -ascii
a = eigenvectors(:,4);
save A10.dat a -ascii
a = eigenvectors(:,5);
save A11.dat a -ascii

save A12.dat eigvals -ascii


% Problem 3

opt_eps = zeros(5,1);
f_vals = zeros(5,1);
xp=-4:.1:4;
K=1;

options = optimoptions('fminunc','Display','iter');

eps_start=[0,2,4,6,8]; % beginning value for epsilon

for modes=1:5 % begin mode loop
    
    eps0=eps_start(modes); % initial value of the eigenvalue epsilon
    f = @(eps)optfunction(eps,xp,K);
    eps = fminunc(f, eps0, options); % fminunc should return epsilon
    opt_eps(modes,1)=eps; % save the eigenvalue
    
end

iterations = [2; 3; 2; 2; 2];

% output answers
save A13.dat opt_eps -ascii
save A14.dat iterations -ascii
