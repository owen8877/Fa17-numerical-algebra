function [H, P] = FrancisQR(H0)
    H = H0;
    n = size(H0, 1);
    m = n-1;
    s = H(m, m) + H(n, n);
    t = H(m, m)*H(n, n) - H(m, n)*H(n, m);
    x = H(1, 1)*H(1, 1) + H(1, 2)*H(2, 1) - s*H(1, 1) + t;
    y = H(2, 1)*(H(1, 1) + H(2, 2) - s);
    z = H(2, 1)*H(3, 2);
    P = eye(n);
    for k = 0:n-3
        [v, b] = householder([x; y; z]);
        Pk = blkdiag(eye(k), eye(3) - b*(v*v'), eye(n-k-3));
        P = P * Pk;
        q = max(1, k);
        H(k+1:k+3, q:n) = H(k+1:k+3, q:n) - (v*v')*H(k+1:k+3, q:n)*b;
        r = min(k+4, n);
        H(1:r, k+1:k+3) = H(1:r, k+1:k+3) - H(1:r, k+1:k+3)*(v*v')*b;
        x = H(k+2, k+1);
        y = H(k+3, k+1);
        if k < n-3
            z = H(k+4, k+1);
        end
    end
    [v, b] = householder([x; y]);
    Pk = blkdiag(eye(n-2), eye(2) - b*(v*v'));
	P = P * Pk;
    H(n-1:n, n-2:n) = H(n-1:n, n-2:n) - (v*v')*H(n-1:n, n-2:n)*b;
    H(1:n, n-1:n) = H(1:n, n-1:n) - H(1:n, n-1:n)*(v*v')*b;
    H = triu(H, -1);
end

