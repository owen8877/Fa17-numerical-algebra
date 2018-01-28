function x = lse(A, b)
    [H, d, R] = QRdecomp(A);
    [m, n] = size(A);
    
    % We first compute Q
    Q = eye(m);
    for i = n:-1:1
        v = H(:, i);
        Q = Q - d(i) * v * (v' * Q);
    end
    
    Q1 = Q(:, 1:n);
    c1 = Q1' * b;
    x = solveUEqn(R, c1);
end

