%% Part 1
% Use Composite Trapezoidal and Composite Simpson's Rules to solve Fredholm
% Integral Equation of the second kind

fprintf("Part 1\n");

%% Section A

fprintf("Section A\n");

a = 0;
b = 1;
f = @(x) x.^2;
K = @(x, t) exp(abs(x-t));
m = 4;

h = (b-a)/m;
x = a + h*(0:m);
x = x';

A = [ h/2*K(x(1:end),x(1)), h*K(x(1:end),x(2:end-1)'), h/2*K(x(1:end),x(end)) ];

A = A - eye(size(A,1));

u = A\-f(x)

%% Section B

fprintf("Section B\n");

a = 0;
b = 1;
f = @(x) x.^2;
K = @(x, t) exp(abs(x-t));
m = 4;

h = (b-a)/m;
x = a + h*(0:m);
x = x';

A = [ h/3*K(x(1:end),x(1)), h/3*K(x(1:end),x(2:end-1)'), h/3*K(x(1:end),x(end)) ];

for i = 2:size(A,1)-1
    if mod(i, 2) == 0
        A(:,i) = 2*A(:,i);
    else
        A(:,i) = 4*A(:,i);
    end
end
A = A - eye(size(A,1));

u = A\-f(x)

%% Section C

fprintf("Section C\n");

a = 0;
b = 1;
f = @(x) x.^2;
K = @(x, t) exp(abs(x-t));
m = 10;

h = (b-a)/m;
x = a + h*(0:m);
x = x';

A = [ h/2*K(x(1:end),x(1)), h*K(x(1:end),x(2:end-1)'), h/2*K(x(1:end),x(end)) ];

A = A - eye(size(A,1));

u = A\-f(x)