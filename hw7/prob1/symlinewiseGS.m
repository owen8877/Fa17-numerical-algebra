function [x, c] = symlinewiseGS(x0, Mblk, Qblk, f, itrLimit, c)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    if ~isfield(c, 'Ml')
        c.Mblk = Mblk;
        c.Qblk = Qblk;
        c.Qblkt = Qblk';
        c.Mblk = Mblk;
        c.nQblk = -Qblk;
        c.nQblkt = c.nQblk';
    end
    
    r0 = f - leftMulMQ(Mblk, Qblk, x0);
    x = x0;
    
    g = solveLEqn_LK(c.Mblk, c.Qblk, f);
    h = solveUEqn_LK(c.Mblk, c.Qblkt, f);
    
    itr = 0;
    while true
        x = solveLEqn_LK(c.Mblk, c.Qblk, leftMulU(c.nQblk, x)) + g;
        x = solveUEqn_LK(c.Mblk, c.Qblkt, leftMulL(c.nQblkt, x)) + h;
        r = f - leftMulMQ(c.Mblk, c.Qblk, x);
        
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

function x = solveUEqn_LK(U, K, y)
    n = size(U, 1);
    x = zeros(n^2, 1);
    
    % x(1+(n-1)*n:n*n) = solveUEqn(U, y(1+(n-1)*n:n*n));
    x(1+(n-1)*n:n*n) = U \ y(1+(n-1)*n:n*n);
    for j = n-1:-1:1
        % x(1+(j-1)*n:j*n) = solveUEqn(U, y(1+(j-1)*n:j*n) - K * x(1+j*n:(j+1)*n));
        x(1+(j-1)*n:j*n) = U \ (y(1+(j-1)*n:j*n) - K * x(1+j*n:(j+1)*n));
    end
end

function b = leftMulL(K, x)
    % the main block diag line is zero
    n = size(K, 1);
    X = reshape(x, n, n);
    KX = K * X;
    b = [zeros(n, 1); reshape(KX(:, 1:n-1), n*(n-1), 1)];
end

function b = leftMulU(K, x)
    % the main block diag line is zero
    n = size(K, 1);
    X = reshape(x, n, n);
    KX = K * X;
    b = [reshape(KX(:, 2:end), n*(n-1), 1); zeros(n, 1)];
end