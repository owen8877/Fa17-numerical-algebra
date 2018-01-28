clc
clear

n = 84;
A = initA(n);
b = initb(n);

xPrecise = A \ b;
xNormalGaussian = normalGaussian(A, b);
xGaussianWithPrimaryColumn = gaussianWithPrimaryColumn(A, b);

figure
s1 = subplot(2, 1, 1);
hold(s1, 'on')
scatter(1:n, abs(xPrecise - xNormalGaussian), 'bo')
scatter(1:n, abs(xPrecise - xGaussianWithPrimaryColumn), 'ro')
title('Absolute Error of GaussElimination with and without Primary Column')
legend('without', 'with')

subplot(2, 1, 2);
scatter(1:n, xPrecise - xGaussianWithPrimaryColumn, 'ro')
title('Error of GaussElimination with Primary Column')

fprintf('Gaussian Method Error\n\t(normal)\t%.3e/%.3e\n\t(improved)\t%.3e/%.3e\nNote: COND(A)=%e\n', ...
    norm(xPrecise - xNormalGaussian), norm(b - A * xNormalGaussian), ...
    norm(xPrecise - xGaussianWithPrimaryColumn), norm(b - A * xGaussianWithPrimaryColumn), ...
    cond(A))

function A = initA(n)
    A = zeros(n, n);
    A(1, 1) = 6;
    for i = 1:n-1
        A(i+1, i+1) = 6;
        A(i, i+1) = 1;
        A(i+1, i) = 8;
    end
end

function b = initb(n)
    b = zeros(n, 1);
    b(1) = 7;
    b(2:n-1) = 15;
    b(n) = 14;
end

function [L, U] = LUDecomposition(A)
    U = A; L = U * 0;
    n = size(U, 1);
    for k = 1:(n-1)
        U(k+1:n, k) = U(k+1:n, k) / U(k, k);
        U(k+1:n, k+1:n) = U(k+1:n, k+1:n) - U(k+1:n, k) * U(k, k+1:n);
    end
    
    for i = 1:n
        L(i, i) = 1;
        for j = 1:i-1
            L(i, j) = U(i, j);
            U(i, j) = 0;
        end
    end
end

function x = normalGaussian(A, b)
    [L, U] = LUDecomposition(A);
    y = solveLEqn(L, b);
    x = solveUEqn(U, y);
end

function x = gaussianWithPrimaryColumn(A, b)
    [P, L, U] = LUDecompositionWithColumnPrimary(A);
    y = solveLEqn(L, suffleB(P, b));
    x = solveUEqn(U, y);
end

function y = solveLEqn(L, b)
    n = size(b, 1);
    y = b;
    for j = 1:n-1
        y(j) = y(j) / L(j, j);
        y(j+1:n) = y(j+1:n) - y(j) * L(j+1:n, j);
    end
    y(n) = y(n) / L(n, n);
end

function x = solveUEqn(U, y)
    n = size(y, 1);
    x = y;
    for j = n:-1:2
        x(j) = x(j) / U(j, j);
        x(1:j-1) = x(1:j-1) - x(j) * U(1:j-1, j);
    end
    x(1) = x(1) / U(1, 1);
end

function [P, L, U] = LUDecompositionWithColumnPrimary(A)
    U = A; L = U * 0;
    n = size(U, 1);
    for k = 1:(n-1)
        [~, primaryColumn] = max(abs(U(k:n, k)));
        primaryColumn = primaryColumn + k - 1;
        P(k) = primaryColumn;
        
        U([primaryColumn, k], :) = U([k, primaryColumn], :);
        
        U(k+1:n, k) = U(k+1:n, k) / U(k, k);
        U(k+1:n, k+1:n) = U(k+1:n, k+1:n) - U(k+1:n, k) * U(k, k+1:n);
    end
    
    for i = 1:n
        L(i, i) = 1;
        for j = 1:i-1
            L(i, j) = U(i, j);
            U(i, j) = 0;
        end
    end
end

function bt = suffleB(P, b)
    bt = b;
    for i = 1:size(P, 2)
        bt([P(i), i]) = bt([i, P(i)]);
    end
end