clear

% define the function
f = @(x) x*sin(3*x)-exp(x);

% define the derivative 
df = @(x) sin(3*x)+x*3*cos(3*x)-exp(x);

% instantiate an array to hold x values
x = zeros(1,1000);

% set the first element of the array
x(1) = -1.6;

% set an error tolerance
tol = 10^(-6);

% newb matlab initiate i to start a counter :(
i=1;

% iterate the newton-raphson method
while abs(f(x(i))) > tol
    x(i+1) = x(i) - f(x(i))/df(x(i));
    i = i+1;
end

% remove all the extra array elements
x = x(1:i);

% output to A1.dat file
save A1.dat x -ascii

% QUESTION 1 PART 2

% define endpoints
a=-.7; b=-.4;

% define a midpoint
c = (a+b)/2;

z = zeros(1,100);
i = 1;

while abs(f(c)) > tol
    z(i) = c;
    if sign(f(c)) == sign(f(b))
        b = c;
        c = (b+a)/2;
    else 
        a = c;
        c = (b+a)/2;
    end
    i=i+1;
end
z(i)=c;
z = z(1:i);

% output to A1.dat file
save A2.dat z -ascii

% output to A3.dat file
w = [length(x) length(z)];
save A3.dat w -ascii

A = [1 2;-1 1];
B = [2 0;0 2];
C = [2 0 -3; 0 0 -1];
D = [1 2; 2 3; -1 0];
x = [1;0]; y = [0;1]; z=[1;2;-1];

% QUESTION 2 PART A
a4 = A + B;
save A4.dat a4 -ascii

% QUESTION 2 PART B
a5 = 3*x - 4*y;
save A5.dat a5 -ascii

% QUESTION 2 PART C
a6 = A*x;
save A6.dat a6 -ascii

% QUESTION 2 PART D
a7 = B*(x-y);
save A7.dat a7 -ascii

% QUESTION 2 PART E
a8 = D*x;
save A8.dat a8 -ascii

% QUESTION 2 PART F
a9 = D*y+z;
save A9.dat a9 -ascii

% QUESTION 2 PART G
a10 = A*B;
save A10.dat a10 -ascii

% QUESTION 2 PART H
a11 = B*C;
save A11.dat a11 -ascii

% QUESTION 2 PART I
a12 = C*D;
save A12.dat a12 -ascii