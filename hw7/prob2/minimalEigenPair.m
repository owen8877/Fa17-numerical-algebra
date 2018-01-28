function [mu, nu] = minimalEigenPair(A, M)
    n = size(A, 1);
%     mu0 = max(abs(diag(A)) - sum(abs(A), 2));
    mu0 = -min(sum(abs(A), 2));
    [mu, x] = inverseMethod(A, 5, mu0);
    nu = (A - sparse(1:n, 1:n, mu*(ones(1, n)))) \ x;
    nu = nu / norm(nu);
    mu = (nu'*A*nu) / (nu'*nu);
end

function [l, x] = powerMethod(A, n)
    x = rand(size(A, 1), 1);
    for j = 1:n
        x = A * x;
        x = x / norm(x);
    end
    l = (x'*A*x) / (x'*x);
end

function [l, x] = inverseMethod(A, itr, mu0)
    n = size(A, 1);
    x = rand(n, 1);
    for j = 1:itr
        x = (A - mu0*eye(n)) \ x;
        x = x / norm(x);
    end
    l = (x'*A*x) / (x'*x);
end
