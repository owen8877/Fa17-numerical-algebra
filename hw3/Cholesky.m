function x = Cholesky(A, b, normal)
    if normal
        x = solve(A, b);
    else
        x = solveImproved(A, b);
    end
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

function L = CholeskyDecomposition(A)
    L = A;
    n = size(A, 1);
    for k = 1:n
        L(k, k) = sqrt(L(k, k));
        L(k+1:n, k) = L(k+1:n, k) / L(k, k);
        for j = k+1:n
            L(j:n, j) = L(j:n, j) - L(j:n, k) * L(j, k);
        end
    end
    
    for i = 1:n
        for j = 1:i-1
            L(j, i) = 0;
        end
    end
end

function [L, D] = ImprovedCholeskyDecomposition(A)
    L = A;
    n = size(A, 1);
    D = zeros(n, n);
    for j = 1:n
        v = zeros(j-1, 1);
        for i = 1:j-1
            v(i) = L(j, i) * L(i, i);
        end
        
        L(j, j) = L(j, j) - L(j, 1:j-1) * v(1:j-1);
        L(j+1:n, j) = (L(j+1:n, j) - L(j+1:n, 1:j-1) * v(1:j-1)) / L(j, j);
    end
    
    for i = 1:n
        D(i, i) = L(i, i);
        for j = 1:i-1
            L(j, i) = 0;
        end
        L(i, i) = 1;
    end
end

function x = solve(A, b)
    L = CholeskyDecomposition(A);
    y = solveLEqn(L, b);
    x = solveUEqn(L', y);
end

function x = solveImproved(A, b)
    [L, D] = ImprovedCholeskyDecomposition(A);
    y = solveLEqn(L, b);
    x = solveUEqn(D * L', y);
end