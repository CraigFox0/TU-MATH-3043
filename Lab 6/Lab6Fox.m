%% Part 1
% Approximates using Composite Trapezoidal and Composite Simpson's Rule

fprintf("Part 1: Approximates using Composite Trapezoidal and Composite Simpson's Rule\n");

a = 0;
b = 2;
n = 10;
h = (b-a)/n;
f = @(x) x^2*log(x^2+1);

syms i;
compositeTrapezoidalApproximation = h*(f(a) + f(b) + 2 * symsum(f(a + (i*h)), i, 1, n-1))/2;
fprintf("Composite Trapezoidal Approximation: %.5f\n", compositeTrapezoidalApproximation);

syms j;
compositeSimpsonApproximation = h*(f(a) + f(b) + 2 * symsum(f(a + 2*j*h), j, 1, (n/2)-1) + 4 * symsum(f(a+(2*j-1)*h), j, 1, n/2))/3;
fprintf("Composite Simpson Approximation: %.5f\n", compositeSimpsonApproximation);

%% Part 2
% Approximates the length of the Gateway Arch

fprintf("Part 2: Approximates the length of the Gateway Arch\n");

a = -91.2;
b = 91.2;
n = 8;
h = (b-a)/n;
diffF = @(x) .03291765*-20.96*sinh(.03291765*x); % Derivate of equation
f = @(x) (1+(diffF(x)).^2).^0.5;

syms j;
approximation = h*(f(a) + f(b) + 2 * symsum(f(a + 2*j*h), j, 1, (n/2)-1) + 4 * symsum(f(a+(2*j-1)*h), j, 1, n/2))/3;
fprintf("Length of Gateway Arch Central Curve using Composite Simpson Approximation: %.0f\n", approximation);

%% Part 3
% Approximate T using the Composite Trapezoidal rule

fprintf("Part 3: Approximate T using the Composite Trapezoidal rule\n");

a = .308;
b = .478;
theta = .7051;
n = 10; %There is 11 points so n = 11 - 1
h = (b-a)/n;
r = [.308, .325, .342, .359, .376, .393, .410, .427, .444, .461, .478];
T = [640, 794, 885, 943, 1034, 1064, 1114, 1152, 1204, 1222, 1239];

A = h*(T(1)*r(1)*theta + T(11)*r(11)*theta + 2 * sum(T(2:n-1).*r(2:n-1)*theta))/2;
B = h*(r(1)*theta + r(11)*theta + 2 * sum(r(2:n-1)*theta))/2;
fprintf("T Approximation: %.4f\n", A/B);