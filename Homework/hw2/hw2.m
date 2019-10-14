% HW2
% Question 1

%Initialization
clear all;
close all;

% set a tolerance for shooting
tol = 10^(-4);

% set eigenfunction colors
col = ['r','b','g','c','m'];


n0=100; % defining parameter n0
A=1; % define initial slope at x=-4
x0=[0 A]; % initial conditions: x1(-4)=0, x1'(-4)=A
xp=[-4 4]; % span of the computational domain

eps_start=n0; % beginning value for epsilon
for modes 1:5 % begin mode loop
    eps=eps_start; % initial value of the eigenvalue epsilon
    deps=n0/100; % default step size in epsilon
    for j=1:1000 % begin convergence loop for beta
        [t,y]=ode45('shoot2',xp,x0,[],n0, eps); % solve ODEs
        if abs(y(end,1)-0) < tol %check for convergence
            disp(eps) % print out the eigenvalue
            break % leave convergence loop
        end
        if (-1)^(modes+1)*y(end,1)>0 % check if epsilon should be higher or lower
            eps=eps-deps; % set epsilon lower by delta epsilon
        else
            eps=eps+deps; % set epsilon higher by delta epsilon
            deps=deps/2; % set the change in epsilon to 1/2 previous change
        end
    end % end converge loop
eps_start=eps-.1; % after finding eigenvalue pick new starting value for next mode
norm=trapz(t,y(:,1).*y(:,1)); % calculate the normalization 
plot(t,y(:,1)/sqrt(norm),col(modes)); hold on % plot each mode


L = 4;
xspan = -L:.1:L;
