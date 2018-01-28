clc
clear

n = 84;
A1 = A1gen(n);
b1 = b1gen(n);
A2 = A2gen(n);
b2 = b2gen(n);

n = 40;
[A3, b3] = Ab3gen(n);

xP1 = A1 \ b1;
xP2 = A2 \ b2;
xP3 = A3 \ b3;

methods = {...
    {'normal gauss', @(A, b) GaussElimination(A, b, true)}, ...
    {'improv gauss', @(A, b) GaussElimination(A, b, false)}, ...
    {'normal chole', @(A, b) Cholesky(A, b, true)}, ...
    {'improv chole', @(A, b) Cholesky(A, b, false)}, ...
    {'QR', @leSolverUsingQR} ...
    };
errorM = zeros(3, 5);

for i = 1:5
    method = methods{i};
    x1 = method{2}(A1, b1);
    errorM(1, i) = norm(x1 - xP1);
    x2 = method{2}(A2, b2);
    errorM(2, i) = norm(x2 - xP2);
    x3 = method{2}(A3, b3);
    errorM(3, i) = norm(x3 - xP3);
end

disp(errorM)

%%
function A = A1gen(n)
    A = zeros(n, n);
    A(1, 1) = 6;
    for i = 1:n-1
        A(i+1, i+1) = 6;
        A(i, i+1) = 1;
        A(i+1, i) = 8;
    end
end

function b = b1gen(n)
    b = zeros(n, 1);
    b(1) = 7;
    b(2:n-1) = 15;
    b(n) = 14;
end

function A = A2gen(n)
    A = zeros(n, n);
    A(1, 1) = 10;
    for i = 1:n-1
        A(i+1, i+1) = 10;
        A(i, i+1) = 1;
        A(i+1, i) = 1;
    end
end

function b = b2gen(n)
    b = rand(n, 1);
    b = b / norm(b);
end

function [A, b] = Ab3gen(n)
   A = zeros(n, n);
   for i = 1:n
       for j = 1:n
           A(i, j) = 1 / (i+j-1);
       end
   end
   b = A * ones(n, 1);
end