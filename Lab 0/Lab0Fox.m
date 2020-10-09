%% Part 2: Letter Grade
% Prints your letter grade based on numerical grade

grade = 94; %set numerical grade here
if grade >= 93
    fprintf("Letter Grade is A\n")
elseif grade >= 90
    fprintf("Letter Grade is A-\n")
elseif grade >= 87
    fprintf("Letter Grade is B+\n")
elseif grade >= 83
    fprintf("Letter Grade is B\n")
elseif grade >= 80
    fprintf("Letter Grade is B-\n")
elseif grade >= 77
    fprintf("Letter Grade is C+\n")
elseif grade >= 73
    fprintf("Letter Grade is C\n")
elseif grade >= 70
    fprintf("Letter Grade is C-\n")
elseif grade >= 67
    fprintf("Letter Grade is D+\n")
elseif grade >= 63
    fprintf("Letter Grade is D\n")
elseif grade >= 60
    fprintf("Letter Grade is D-\n")
else
    fprintf("Letter Grade is F\n")
end
%% Part 3: Square Roots
% Prints a table of first ten non-zero, positive integers and their square roots

fprintf('Number, Squareroot\n')
for num = 1:10
    fprintf('%d,\t%f \n', num, num^0.5)
end
%% Part 4: Partial Sum
% Calculates the partial sum for (2/3)^k for n = 5, 10, 25, 50, 100

total = 0;
for n = 1:100
    total = total + (2/3)^n;
    if n == 5 | n == 10 | n == 25 | n == 50 | n == 100
        fprintf('The partial sum for %d is %f\n', n, total)
    end
end
%% Part 5: Divder
% Divides a number and outputs result until it is less than 1

num = 843; %Set number to be divided here
while num > 1
    num = num / 2;
    fprintf('%f\n', num)
end