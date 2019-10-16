% main program u''(x) = -sin(u)
close all; clear all; clc;
% create params
params.x0 = 0;
params.xf = round(2*pi, 3);  % being smart 
params.ffunc = @(theta) -sin(theta);
params.ffuncDer = @(theta) -cos(theta);
%params.ffunc = @(theta) -2*sin(theta) + exp(-0.5*theta);
%params.ffuncDer = @(theta) -2*cos(theta) - 0.5*exp(-0.5*theta);
params.h = 1e-3; 

params.alpha = 0.7; 
params.beta = 0.7; 

%%
N = (params.xf-params.x0)/params.h; % number of steps 
assert(N-round(N) == 0); % make sure number of steps is even. 
xs = linspace(params.x0, params.xf, N); 

% initialization function and its derivative (for bvp4c)
theta_in = @(x) 0.7*cos(x) - 0.5*sin(x); 
theta_in_der = @(x) -0.7*sin(x) - 0.5*cos(x);
params.infunc = theta_in; 
params.infuncDer = theta_in_der; 

% initialize
theta = theta_in(xs); 
theta = theta'; % has to be a column vector 

[G,J] = buildSystemNL(theta, params);
% measure initial error 
err = norm(G);
iter = 0; 
thetas = []; % lets store these
while(err > 1e-7) 
    % take Newton step
    iter = iter+1;
    deltaN = -J\G;    
    deltaGN = -(J'*J)\(J'*G);
%    deltaGN = -(J'*J+0.1*speye(size(J'*J)))\(J'*G);
    delta = deltaGN;
    eold = norm(G);
    enew = eold + 1; 
    step = 2;
    while(enew > 0.9*eold)
        step = step/2;
        thetac = theta + step*delta;
        [G,J] = buildSystemNL(thetac, params);
        enew = norm(G);
        if step <1e-10
            error('step too small, failed to solve')
        end
    end
    theta = thetac;
    thetas = [thetas, theta]; % add a column
%    [G,J] = buildSystemNL(theta, params);
    err = enew;
    fprintf('iter: %d, error: %7.2d, step: %7.3e\n', iter, err, step);
    if(iter > 100)
        break;
    end
end

plot(theta, 'Linewidth', 2)
hold on
plot(thetas(:,1), 'Linewidth', 1);
plot(thetas(:,2), 'Linewidth', 1);
plot(thetas(:,3), 'Linewidth', 1);

legend('final sol', '1 iter', '2 iters','3 iters' );


%% Checking our work with bvp4c 
sol = bvp4cRunner(params); 

%%
figure(2)
plot(theta, 'Linewidth', 2)
hold on
plot(sol.y(1,:), '--', 'Linewidth', 2);
legend('our sol', 'bvp4c sol');

%% check bvp4c's solution
[G,J] = buildSystemNL(sol.y(1,:)', params);
norm(G(2:end-1))