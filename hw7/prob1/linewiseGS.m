function x = linewiseGS(x0, Mblk, Qblk, f, itrLimit)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    r0 = f - leftMulMQ(Mblk, Qblk, x0);
    x = x0;
    
    Ml = tril(Mblk);
    nMu = -triu(Mblk, 1);
    nQblk = -Qblk;
    g = solveLEqn_LK(Ml, Qblk, f);
    
    itr = 0;
    while true
        x = solveLEqn_LK(Ml, Qblk, leftMulU(nMu, nQblk, x)) + g;
        r = f - leftMulMQ(Mblk, Qblk, x);
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
        
        if (norm(r) < norm(r0) * 1e-6) || itr >= itrLimit
            break
        end
    end
end

function x = solveLEqn_LK(L, K, y)
    n = size(L, 1);
    x = zeros(n^2, 1);
    
    % x(1:n) = solveLEqn(L, y(1:n));
    x(1:n) = L \ y(1:n);
    for j = 2:n
        % x(1+(j-1)*n:j*n) = solveLEqn(L, y(1+(j-1)*n:j*n) - K * x(1+(j-2)*n:(j-1)*n));
        x(1+(j-1)*n:j*n) = L \ (y(1+(j-1)*n:j*n) - K * x(1+(j-2)*n:(j-1)*n));
    end
end

function b = leftMulU(U, K, x)
    n = size(U, 1);
    X = reshape(x, n, n);
    UX = U * X;
    KX = K * X;
    UX(:, 1:n-1) = UX(:, 1:n-1) + KX(:, 2:end);
    b = reshape(UX, n^2, 1);
end