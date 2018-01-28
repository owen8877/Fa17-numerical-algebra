function [x, history] = naiveSOR(A, b)
    D = diag(diag(A));
    L = -A; U = -A;
    n = size(A, 1);
    for i = 1:n
        for j = 1:i
            L(j, i) = 0;
            U(i, j) = 0;
        end
    end
    omega = 1.3;
    M = (D-omega*L) \ ((1-omega)*D + omega*U);
    g = omega * ((D-omega*L) \ b);
    [x, history] = iterFramework(M, g, zeros(n, 1));
end
