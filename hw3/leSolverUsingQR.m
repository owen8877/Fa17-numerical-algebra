function x = leSolverUsingQR(A, b)
    [H, d, R] = QRdecomp(A);
    n = size(A, 1);
    
    for i = 1:n
        v = H(:, i);
        % v(i) = 1;
        b = b - d(i) * v * (v' * b);
    end
    
    x = solveUEqn(R, b);
end

