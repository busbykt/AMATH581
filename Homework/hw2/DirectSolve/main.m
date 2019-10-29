% main program to solve u''(x) = f(x)
close all; clear all; clc;
% create params
params.x0 = 0;
params.xf = 1; 
params.ffunc = @(x) sin(10*x); 
params.h = .01; 
params.alpha = 0; 
params.beta = 0; 

[A,F] = buildSystem(params);

% compute solution
u = A\F; 

figure()
plot(u);

%figure();
%uxx = A*u; 
%plot(uxx);



