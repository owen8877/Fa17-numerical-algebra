clc
clear

load data.mat

A = [ones(size(data, 1), 1) data];
b = price;

x = lse(A, b);
xP = (A'*A) \ (A'*b);
fprintf('LS Result:\n');
disp(x);
fprintf('Error\t%e\n', norm(x - xP));