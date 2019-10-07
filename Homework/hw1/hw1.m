clear

% QUESTION 1 PART A

% Define the ODE
dy = @(y,t) -3*y(t)*sin(t);
% Define the analytical so  lution
y_true = @(t) pi()*exp(3*(cos(t)-1))/sqrt(2);

% Forward Euler Implementation
yFE = zeros(7,10000);
Y = zeros(7,10000);
err = zeros(1,7);
yFE(1:7) = pi/sqrt(2);
Y(1:7) = y_true(0);
del_t = zeros(1,7);

% Heun Method Implementation
yH = zeros(7, 10000);
errh = zeros(1,7);
yH(1:7) = pi/sqrt(2);
fH = @(t,y) -3*y*sin(t);

i = 1;
while i+1 <= 8
    t_vals = 0:2^-(i+1):5;
    del_t(i) = 2^-(i+1);
    j = 1;
    for t_val = t_vals
        yFE(i,j+1) = yFE(i,j) + del_t(i)*-3*yFE(i,j)*sin(t_val);
        yH(i,j+1) = yH(i,j) + del_t(i)/2*(-3*yH(i,j)*sin(t_val) + -3*(yH(i,j)+del_t(i)*-3*yH(i,j)*sin(t_val))*sin(t_val+del_t(i)));
        %yH(i,j+1) = yH(i,j) + del_t(i)/2*(-3*yH(i,j)*sin(t_vals(j)) + f(t_vals(j+1),yH(i,j)+del_t(i)*-3*yH(i,j)*sin(t_vals(j))));
        Y(i,j) = y_true(t_val);
        j = j + 1;
    end
    j = length(t_vals);
    % calculate error for the given t_val
    err(i) = mean(abs(Y(i,1:j) - yFE(i,1:j)));
    errh(i) = mean(abs(Y(i,1:j) - yH(i,1:j)));
    i = i + 1; 
    
end

% save the output to a column vector
sol = transpose(yFE(7,1:j));
solh = transpose(yH(7,1:j));

% output 2e-8 step yFE & yH to A1.dat & A4.dat files
save A1.dat sol -ascii
save A4.dat solh -ascii

% output error values for all steps to A2.dat
save A2.dat err -ascii
save A5.dat errh -ascii

scatter(log(del_t), log(err))

p = polyfit(log(del_t), log(err),1);
slope = p(1);   
p_data = polyval(p, log(del_t));
hold on 
plot(log(del_t), p_data)

% output linear polyfit slope to A3.dat
save A3.dat slope -ascii

% QUESTION 1 PART B
hold on
scatter(log(del_t), log(errh))
title('Log Error Euler & Heun over Log Delta T')     
xlabel('Log Delta T')
ylabel('Log Error')

p = polyfit(log(del_t), log(errh),1);
slopeh = p(1);   
p_data = polyval(p, log(del_t));
hold on 
plot(log(del_t), p_data)

% output linear polyfit slope to A3.dat
save A6.dat slopeh -ascii

% QUESTION 2

eps = .1;
t_vals = 0:.5:32;
[t,y] = ode45(@(t,y) vanderpol(t,y,eps),t_vals,[sqrt(3); 1]);
hold on
plot(t,y(:,1),'-o',t,y(:,2),'-o')
title("Solution of van der Pol Equation (\epsilon ="+ eps+") with ODE45");
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

eps = 1;
[t,yy] = ode45(@(t,yy) vanderpol(t,yy,eps),t_vals,[sqrt(3); 1]);
hold on
plot(t,yy(:,1),'-o',t,yy(:,2),'-o')
title("Solution of van der Pol Equation (\epsilon ="+ eps+") with ODE45");
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

eps = 20;
[t,yyy] = ode45(@(t,yyy) vanderpol(t,yyy,eps),t_vals,[sqrt(3); 1]);
hold on
plot(t,yyy(:,1),'-o',t,yyy(:,2),'-o')
title("Solution of van der Pol Equation (\epsilon ="+ eps+") with ODE45");
xlabel('Time t');
ylabel('Solution y');
legend('y_1','y_2')

solution = [y(:,1) yy(:,1) yyy(:,1)];

% output linear polyfit slope to A8.dat
save A7.dat solution -ascii


% QUESTION 2 PART B 

eps = 1;
tspan = [0 32];
y0 = [2 ; pi^2];
TOLS = logspace(-4,-10,7);
del_t = zeros(1,length(TOLS));

i = 1;
for TOL = TOLS
    options = odeset('AbsTol',TOL,'RelTol',TOL); 
    [T,Y] = ode45(@(t,y) vanderpol(t,y,eps),tspan,y0,options);
    del_t(i) = mean(diff(T));
    i = i+1;
end 

scatter(log(del_t), log(TOLS))
xlabel('Log Delta t')
ylabel('Log Tolerance')

p = polyfit(log(del_t), log(TOLS), 1);
sol = p(1);

% output linear polyfit slope to A8.dat
save A8.dat sol -ascii

%ODE 23 SOLUTION

del_t = zeros(1,length(TOLS));

i = 1;
for TOL = TOLS
    options = odeset('AbsTol',TOL,'RelTol',TOL); 
    [T,Y] = ode23(@(t,y) vanderpol(t,y,eps),tspan,y0,options);
    del_t(i) = mean(diff(T));
    i = i+1;
end 

scatter(log(del_t), log(TOLS))
xlabel('Log Delta t')
ylabel('Log Tolerance')

p = polyfit(log(del_t), log(TOLS), 1);
sol = p(1);

% output linear polyfit slope to A9.dat
save A9.dat sol -ascii

% ODE 113 SOLUTION

del_t = zeros(1,length(TOLS));

i = 1;
for TOL = TOLS
    options = odeset('AbsTol',TOL,'RelTol',TOL); 
    [T,Y] = ode113(@(t,y) vanderpol(t,y,eps),tspan,y0,options);
    del_t(i) = mean(diff(T));
    i = i+1;
end 

scatter(log(del_t), log(TOLS))
xlabel('Log Delta t')
ylabel('Log Tolerance')

p = polyfit(log(del_t), log(TOLS), 1);
sol = p(1);

% output linear polyfit slope to A10.dat
save A10.dat sol -ascii




