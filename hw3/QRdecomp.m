function [H, d, R] = QRdecomp(A)
    [m, n] = size(A);
    
    d = zeros(n, 1);
    H = zeros(m, n);
    for j = 1:n
        if j < m
            [v, b] = householder(A(j:m, j));
            A(j:m, j:n) = A(j:m, j:n) - b * v * (v' * A(j:m, j:n));
            d(j) = b;
            A(j+1:m, j) = v(2:m-j+1);
            H(j:m, j) = v(1:m-j+1);
        end
    end
    
    R = zeros(m, n);
    for j = 1:n
        for i = 1:min(j, m)
            R(i, j) = A(i, j);
        end
    end
end

