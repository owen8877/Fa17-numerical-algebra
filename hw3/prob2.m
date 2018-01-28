clc
clear

t = [-1 -0.75 -0.5 0 0.25 0.5 0.75]';
A = [t.^2 t ones(7, 1)];
b = [1 0.8125 0.75 1 1.3125 1.75 2.3125]';

x = lse(A, b);
xP = (A'*A) \ (A'*b);
fprintf('LS Result:\n\ta\t%e\n\tb\t%e\n\tc\t%e\n', x(1), x(2), x(3));
fprintf('Error\t%e\n', norm(x - xP));