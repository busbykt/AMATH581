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
y0=[1 4]; % initial conditions: y1(-4)=1, y1'(-4)=4
xp=-4:.1:4; % span of the computational domain

eigenvalues = zeros(1,5);
eigenvectors = zeros(81,5);

eps_start=K; % beginning value for epsilon
for modes=1:5 % begin mode loop
    eps=eps_start; % initial value of the eigenvalue epsilon
    deps=K/100; % default step size in epsilon
    for j=1:1000 % begin convergence loop for beta
        [t,y]=ode45('shoot2',xp,y0,[],K, eps); % solve ODEs
        if abs(y(end,1)-0) < tol %check for convergence
            disp(eps) % print out the eigenvalue
            eigenvalues(1,modes)=eps;
            eigenvectors(:,modes)=y(:,1);
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
end
legend('mode1','mode2','mode3','mode4','mode5');
