%% Defining an anonymous function
% One of the easiest ways to define a function that takes in variables is
% via anonymous functions.
%
% Suppose we want to define a function f(x) = sin(x). We would define an
% anonymous funciton for f in the following way:
 f = @(x) sin(x); 
% The syntax here means: 
% f = function name
% @(x) = input variable(s) are only x 
% sin(x) = defined function of x

% Anonymous functions can be defined using several variables as well. See 
% https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html
% for more details on anonymous functions.

%% Vectorized code & Elementwise operations
% One very convenient aspect about MATLAB is its ability to perform
% elementwise operations on vectors/matrices, which allows for vectorized
% code. 
%
% Note that for vectors/matrices, addition and subtraction of other
% vectors/matrices is already an elementwise operation in theory. However,
% there is no mathematical definition of a vector times vector operation
% that returns a vector, and a matrix multiplication is not an elementwise
% operation. 
% As we will see in the section on plotting, somethimes we want to perform 
% these elementwise operations, so we can use '.*' to perform elementwise
% multiplication, './' to perform elementwise division, and '.^' to perform
% elementwise powers.

% Let's see a few examples:
a = [1;2;3;4]; b = [5;6;7;8]; 
% multiply vectors a and b elementwise, output should be [5;12;21;32]
fprintf('\na.*b = \n');
c = a.*b
% Divide b by a elementwise, output should be [5;3;2.333333;2]
fprintf('\nb./a = \n');
d = b./a
% Square each element of a, output should be [1;4;9;16]
fprintf('\na.^2 = \n');
e = a.^2

% Now, we can use elementwise operations on vectors/matrices to vectorize
% code; i.e. use vectors and vector operations to define a new vector
% without using loops.
% Below is an example of vectorized vs non-vectorized code for finding a
% vector y containing the squares of the first 5 integers:

% Non-vectorized:
y = zeros(5,1);
for int = 1:5
    y(int) = int^2;
end
% Vectorized:
int = 1:5; % define int as a vector of integers
y = int.^2; % square each integer

% Generally, vectorized code is preferable to loop-based code whenever
% possible, as it is cleaner to read and usually runs faster.

%% Plotting basics
% In short, MATLAB's plotter works by the user specifying a vector of x
% coordinates along with a corresponding vector of y coordinates, and then
% MATLAB plots the graph of these coordinates on the (x,y) plane. 
%
% Let's plot the first 5 integers and their squares, using 'int' from the
% previous example as the x values, and 'y' as the y values:
figure % produces a new plotting screen
plot(int,y)

% As we can see, the line looks jagged. This is because MATLAB's plotter
% plots the points on the (x,y) plane that we specify, and then draws lines
% between the points to connect the dots. If we want a higher resolution
% plot, we need to use more coordinates:
x = linspace(1,5,1000); % generates 1000 linearly spaced points between 1 and 5
y = x.^2; 
figure
plot(x,y,'*r') % plots the new plot in red with '*' stars marking the points


% now, let's plot 2 things on the same screen:
f1 = @(x) sin(x); % sin(x)
f2 = @(x) sin(x.^2); % sin(x^2)
x = linspace(0,2*pi,1000); % generates 1000 linearly spaced points between 0 and 2*pi

figure % create new plotting window
plot(x,f1(x)); %plot sin(x) first
hold on % don't clear the plotting window before plotting the next part
plot(x,f2(x),'r'); % plot sin(x.^2) second using a red line
xlabel('x'); ylabel('y'); title('sin(x) vs sin(x^2)'); % label your plot
% give your plot a legend. Legends assign labels in the order of plotting.
legend('sin(x)','sin(x^2)'); 

% See https://www.mathworks.com/help/matlab/ref/plot.html for more details
% on MATLAB's plotter