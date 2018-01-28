function H = Hessenberg(A)
    n = size(A, 1);
    for k = 1:n-2
        [v, b] = householder(A(k+1:n, k));
        A(k+1:n, k:n) = A(k+1:n, k:n) - (v*v')*A(k+1:n, k:n)*b;
        A(1:n, k+1:n) = A(1:n, k+1:n) - A(1:n, k+1:n)*(v*v')*b;
    end
    H = triu(A, -1);
end

