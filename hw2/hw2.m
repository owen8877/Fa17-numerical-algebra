clc
clear

%% estimate \kappa(A) where A is n-order hilbert matrix
kappa = zeros(20, 1);
actual = zeros(20, 1);
for n = 5:20
    A = hilbert(n);
    [L, U] = LUDecomposition(A);
    kappa(n) = approxInfNorm(L, U) * infNorm(A);
    actual(n) = norm(A, Inf) * norm(invhilb(n), Inf);
    fprintf('%d\t%e\t%e\n', n, kappa(n), actual(n));
end

%% second problem
for n = 5:30
    x = rand(n, 1);
    x = x / norm(x);
    A = secondA(n);
    b = A * x;
    
    [P, L, U] = LUDecompositionWithColumnPrimary(A);
    xhat = gaussianWithPrimaryColumn(A, b);
    nu = approxInfNormWithP(P, L, U);
    gamma = max(abs(b - A * xhat));
    beta = max(abs(b));
    mu = infNorm(A);
    rou = nu * mu * gamma / beta;
    fprintf('%d\t%e\t%e\n', n, rou, max(abs(x-xhat))/max(abs(x)));
end

function A = hilbert(n)
    A = zeros(n, n);
    for i = 1:n
        for j = 1:n
            A(i, j) = 1 / (i + j - 1);
        end
    end
end

function A = secondA(n)
    A = eye(n);
    for i = 1:n-1
        for j = 1:i
            A(i+1, j) = -1;
        end
        A(i, n) = 1;
    end
end

function Mnorm = infNorm(A)
    Mnorm = max(sum(abs(A), 2));   
end

function nu = approxInfNorm(L, U)
    n = size(L, 1);
    x = ones(n, 1) / n;
    while true
        w = solveUEqn(L', solveLEqn(U', x));
        v = sign(w);
        z = solveUEqn(U, solveLEqn(L, v));
        if max(abs(z)) <= z' * x
            nu = norm(w, 1);
            break
        else
            x = x * 0;
            x(abs(z) == max(abs(z))) = 1;
        end
    end
end

function nu = approxInfNormWithP(P, L, U)
    n = size(L, 1);
    x = ones(n, 1) / n;
    while true
        w = shuffleB(P, solveUEqn(L', solveLEqn(U', x)));
        v = sign(w);
        z = solveUEqn(U, solveLEqn(L, shuffleB(P, v)));
        if max(abs(z)) <= z' * x
            nu = norm(w, 1);
            break
        else
            x = x * 0;
            x(abs(z) == max(abs(z))) = 1;
        end
    end
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
    y = solveLEqn(L, shuffleB(P, b));
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

function bt = shuffleB(P, b)
    bt = b;
    for i = 1:size(P, 2)
        bt([P(i), i]) = bt([i, P(i)]);
    end
end
