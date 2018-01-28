function y = solveLEqn(L, b)
    n = size(b, 1);
    y = b;
    for j = 1:n-1
        y(j) = y(j) / L(j, j);
        y(j+1:n) = y(j+1:n) - y(j) * L(j+1:n, j);
    end
    y(n) = y(n) / L(n, n);
end
