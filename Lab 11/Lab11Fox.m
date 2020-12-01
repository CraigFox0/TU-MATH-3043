%% Dummy
% Ignore this. It is used to make the outputs allign for the main code

A = [10, -1, 2, 0; -1, 11, -1, 3; 2, -1, 10, -1; 0, 3, -1, 8];
b = [6; 25; -11; 15];
N = 100;
tolerance = 10^-3;
x0 = zeros(size(b, 1), 1);

[jacobiSolution, jacobiIterations] = Jacobi(A, b, N, tolerance, x0);
[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);

%% Part 2
% Use the Jacobi method and the Gauss-Seidel method to solve the linear system

fprintf("Part 2: Use the Jacobi method and the Gauss-Seidel method to solve the linear system\n");

A = [10, -1, 2, 0; -1, 11, -1, 3; 2, -1, 10, -1; 0, 3, -1, 8];
b = [6; 25; -11; 15];
N = 100;
tolerance = 10^-3;
x0 = zeros(size(b, 1), 1);

fprintf("Part A: Jacobi method\n");
[jacobiSolution, jacobiIterations] = Jacobi(A, b, N, tolerance, x0);
fprintf('Jacobi Method Solution\n');
fprintf('%.4f\n', jacobiSolution);
fprintf('Jacobi Iterations: %.0f\n', jacobiIterations);
fprintf("Part B: Gauss-Seidel method\n");
[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution\n');
fprintf('%.4f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

%% Part 3
% Use the Jacobi method and the Gauss-Seidel method to solve the linear system

fprintf("Part 3: Use the Jacobi method and the Gauss-Seidel method to solve the linear system\n");

A = [4, -1, 0, 0, 0, 0; -1, 4, -1, 0, 0, 0; 0, -1, 4, 0, 0, 0; 0, 0, 0, 4, -1, 0; 0, 0, 0, -1, 4, -1; 0, 0, 0, 0, -1, 4];
b = [0; 5; 0; 6; -2; 6];
N = 100;
tolerance = 10^-4;
x0 = zeros(size(b, 1), 1);

fprintf("Part A: Jacobi method\n");
[jacobiSolution, jacobiIterations] = Jacobi(A, b, N, tolerance, x0);
fprintf('Jacobi Method Solution\n');
fprintf('%.5f\n', jacobiSolution);
fprintf('Jacobi Iterations: %.0f\n', jacobiIterations);
fprintf("Part B: Gauss-Seidel method\n");
[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution\n');
fprintf('%.5f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

%% Part 4
% Use the best method to solve the linear system

fprintf("Part 4: Use the best method to solve the linear system\n");

A = [2, -1, 1; 2, 2, 2; -1, -1, 2];
b = [-1; 4; -5];
N = 100;
tolerance = 10^-3;
x0 = zeros(size(b, 1), 1);

fprintf("Part A: Calculate Spectral Radius\n");
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
Tj = inv(D)*(L+U);
Tg = inv(D-L)*U;
pTj = max(abs(eig(Tj)));
fprintf('p(Tj) is %.4f\n', pTj);
pTg = max(abs(eig(Tg)));
fprintf('p(Tg) is %.4f\n', pTg);

fprintf("Since p(Tg) is less than 1 and p(Tj) is not, it is best to use the Gauss-Seidel method.\n");

fprintf("Part B: Solve\n");
[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution\n');
fprintf('%.4f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

%% Part 5
% Use the best method to solve the linear system

fprintf("Part 5: Use the best method to solve the linear system\n");

A = [1, 2, -2; 1, 1, 1; 2, 2, 1];
b = [7; 2; 5];
N = 100;
tolerance = 10^-3;
x0 = zeros(size(b, 1), 1);

fprintf("Part A: Calculate Spectral Radius\n");
D = diag(diag(A));
L = -tril(A,-1);
U = -triu(A,1);
Tj = inv(D)*(L+U);
Tg = inv(D-L)*U;
pTj = max(abs(eig(Tj)));
fprintf('p(Tj) is %.4f\n', pTj);
pTg = max(abs(eig(Tg)));
fprintf('p(Tg) is %.4f\n', pTg);

fprintf("Since p(Tj) is less than 1 and p(Tg) is not, it is best to use the Jacobi method.\n");

fprintf("Part B: Solve\n");
fprintf("Part A: Jacobi method\n");
[jacobiSolution, jacobiIterations] = Jacobi(A, b, N, tolerance, x0);
fprintf('Jacobi Method Solution\n');
fprintf('%.4f\n', jacobiSolution);
fprintf('Jacobi Iterations: %.0f\n', jacobiIterations);

%% Part 6
% Calculate the probability of reaching the left endpoint before the right

fprintf("Part 6: Calculate the probability of reaching the left endpoint before the right\n");

fprintf("Part A: alpha = 1/2\n");

N = 1000;
tolerance = 10^-10;

n = 10;
A = eye(n) - diag(.5*ones(n-1,1), -1) - diag(.5*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = 1/2;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

n = 50;
A = eye(n) - diag(.5*ones(n-1,1), -1) - diag(.5*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = 1/2;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

n = 100;
A = eye(n) - diag(.5*ones(n-1,1), -1) - diag(.5*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = 1/2;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

fprintf("Part B alpha = 1/3\n");

alpha = 1/3;

n = 10;
A = eye(n) - diag(alpha*ones(n-1,1), -1) - diag((1-alpha)*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = alpha;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

n = 50;
A = eye(n) - diag(alpha*ones(n-1,1), -1) - diag((1-alpha)*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = alpha;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

n = 100;
A = eye(n) - diag(alpha*ones(n-1,1), -1) - diag((1-alpha)*ones(n-1,1), 1);
b = zeros(n, 1);
b(1) = alpha;
x0 = zeros(size(b, 1), 1);

[gaussSeidelSolution, gaussSeidelIterations] = GaussSeidel(A, b, N, tolerance, x0);
fprintf('Gauss-Seidel Method Solution for n = %.0d\n', n);
fprintf('%.12f\n', gaussSeidelSolution);
fprintf('Gauss-Seidel Iterations: %.0f\n', gaussSeidelIterations);

%% Local Functions

function [solution, iterations] = Jacobi(A, b, N, tolerance, x0)
answerFound = false;
n = size(b,1);
x = zeros(n, 1);
k = 1;
while k <= N && ~answerFound
    for i = 1:n
        x(i) = (1/A(i,i))*(-A(i,:)*x0(:) + A(i,i)*x0(i) + b(i));
    end
    if norm(x-x0, inf)/norm(x, inf) < tolerance
        iterations = k;
        answerFound = true;
    end
    x0 = x;
    k = k + 1;
end
if answerFound == false
    iterations = N;
end
solution = x;
end

function [solution, iterations] = GaussSeidel(A, b, N, tolerance, x0)
answerFound = false;
n = size(b,1);
x = zeros(n, 1);
k = 1;
while k <= N && ~answerFound    
    for i = 1:n
        x(i) = (1/A(i,i))*(-A(i,1:i-1)*x(1:i-1) - A(i,i+1:end)*x0(i+1:end) + b(i));
    end
    if norm(x-x0, inf)/norm(x, inf) < tolerance
        iterations = k;
        answerFound = true;
    end
    x0 = x;
    k = k + 1;
end
if answerFound == false
    iterations = N;
end
solution = x;
end