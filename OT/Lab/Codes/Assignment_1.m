%% QUESTION 1
clc
clear all
n = 10; % Define the number of elements in the array
natural_numbers = 1:n; % Creates an array of natural numbers from 1 to n
% Display the generated array
disp(natural_numbers);

%% QUESTION 2
clc
clear all
format short
% #1: input parameter
A=[4 1 3; 2 6 7; 3 1 8];
x=min(min(A))
[i,j]=find(A==min(min(A)))
A(i,j)=10

%% QUESTION 3
clc
clear all
format short
% #1: input parameter
A=[4 1 3; 2 6 7; 3 1 8];
% Sort each row in ascending order
sorted_rows = sort(A, 2);
% Sort each column in ascending order
sorted_columns = sort(A, 1);

A= sort(A,2);
A= sort(A,1);

disp("Sorted Rows:");
disp(sorted_rows);

disp("Sorted Columns:");
disp(sorted_columns);

disp("Final sorted matrix");
disp(A);

%% QUESTION 4
clc
clear all
x=0:0.01:12; % values for x for the plot
% the function
y = (12 - 2*x)/4;
% Plotting
plot(x,y,'r.')
xlabel('x \rightarrow')
ylabel('2x+4y=12 \rightarrow')
title('Graph of x vs y')
grid on

%% QUESTION 5
clc
clear all
% The equations of the lines
% eq1 => 2*x + 4*y == 12;
% eq2 => 3*x + 2*y == 12;

% Define the coefficient matrix and the constant vector
A = [2 4; 3 2];
B = [12; 12];
% Solve the system of equations using matrix operations
solution = inv(A)*B;
% Display the intersection point
disp('Intersection Point:');
disp (solution);


