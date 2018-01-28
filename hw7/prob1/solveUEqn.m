function x = solveUEqn(U, y)
    n = size(y, 1);
    x = y;
    for j = n:-1:2
        x(j) = x(j) / U(j, j);
        x(1:j-1) = x(1:j-1) - x(j) * U(1:j-1, j);
    end
    x(1) = x(1) / U(1, 1);
end
