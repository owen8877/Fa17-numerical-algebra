function [x, history] = SOR(A, b)
    D = diag(diag(A));
    L = -A; U = -A;
    n = size(A, 1);
    for i = 1:n
        for j = 1:i
            L(j, i) = 0;
            U(i, j) = 0;
        end
    end
    B = diag(diag(A).^-1) * (L+U);
    % rou = powermethod(B);
    rou = max(abs(eig(B)));
    omega = 2 / (1+sqrt(1-rou^2));
    fprintf('(debug)\tSOR: Best omega is %.4f\n', omega);
    M = (D-omega*L) \ ((1-omega)*D + omega*U);
    g = omega * ((D-omega*L) \ b);
    [x, history] = iterFramework(M, g, zeros(n, 1));
end
