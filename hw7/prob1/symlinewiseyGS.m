function [x, c] = symlinewiseyGS(x0, Mblk, Qblk, f, itrLimit, c)
    % suppose that A has Mblk on the main diag line, with Qblk on the
    % quasi-diag line.
    
    % first construct the whole sparse matrix A
    if ~isfield(c, 'PAP')
        n = size(Mblk, 1);
%         [mi, mj, mv] = find(Mblk);
%         [qi, qj, qv] = find(Qblk);
%         nmblk = numel(mi);
%         nqblk = numel(qi);

%         av = [repmat(mv', 1, n) repmat(qv', 1, 2*(n-1))];
%         ai = [repmat(mi', 1, n) repmat(qi', 1, 2*(n-1))];
%         aj = [repmat(mj', 1, n) repmat(qj', 1, 2*(n-1))];
%         mshift = reshape(repmat(0:n-1, nmblk, 1), 1, n*nmblk);
%         qushift = reshape(repmat(0:n-2, nqblk, 1), 1, (n-1)*nqblk);
%         qlshift = reshape(repmat(1:n-1, nqblk, 1), 1, (n-1)*nqblk);
%         ai = ai + [mshift qushift qlshift] * n;
%         aj = aj + [mshift qlshift qushift] * n;
%         A = sparse(ai, aj, av);

        % then we arrange A in a special way
        indexBase = reshape(1:n^2, n, n);
        d = diag(indexBase);
        [~, ~, u] = find(triu(indexBase, 1)');
        [~, ~, l] = find(tril(indexBase, -1));
%         P = sparse([d; u; l], [d; l; u], ones(1, n^2));

%         PAP = P * A * P;
        
%         D = sparse(n^2, n^2);
%         L = sparse(n^2, n^2);
%         U = sparse(n^2, n^2);

%         for i = 1:n
%             D(1+(i-1)*n:i*n, 1+(i-1)*n:i*n) = PAP(1+(i-1)*n:i*n, 1+(i-1)*n:i*n);
%             L(1+(i-1)*n:i*n, 1:(i-1)*n) = -PAP(1+(i-1)*n:i*n, 1:(i-1)*n);
%             U(1+(i-1)*n:i*n, 1+i*n:end) = -PAP(1+(i-1)*n:i*n, 1+i*n:end);
%         end
        
%         yMblk = PAP(1:n, 1:n);
%         yQblk = PAP(1+n:2*n, 1:n);
        m1 = Mblk(1, 1); m2 = Mblk(1, 2);
        q1 = Qblk(1, 1); q2 = Qblk(1, 2);
        yMblk = createTriDiagSparse(m1, q1, n);
        yQblk = createTriDiagSparse(m2, q2, n);
        
        c.Mblk = yMblk;
        c.Qblk = yQblk;
        c.Qblkt = yQblk';
        c.Mblk = yMblk;
        c.nQblk = -yQblk;
        c.nQblkt = c.nQblk';
%         c.P = P;
        c.dul = [d; u; l];
        c.dlu = [d; l; u];
%         c.D = D;
%         c.L = L;
%         c.U = U;
    end
    
    %Pf = c.P * f;
    Pf(c.dul) = f(c.dlu); Pf = Pf';
    Px(c.dul) = x0(c.dlu); Px = Px';
    r0 = Pf - leftMulMQ(c.Mblk, c.Qblk, Px);
%     r0 = Pf - PAP * Px;
    
    g = solveLEqn_LK(c.Mblk, c.Qblk, Pf);
    h = solveUEqn_LK(c.Mblk, c.Qblkt, Pf);

%     g = (c.D-c.L) \ Pf;
%     h = (c.D-c.U) \ Pf;
    
    itr = 0;
    while true
        Px = solveLEqn_LK(c.Mblk, c.Qblk, leftMulU(c.nQblk, Px)) + g;
        Px = solveUEqn_LK(c.Mblk, c.Qblkt, leftMulL(c.nQblkt, Px)) + h;
        r = Pf - leftMulMQ(c.Mblk, c.Qblk, Px);
%         Px = (c.D-c.L) \ (c.U*Px) + g;
%         Px = (c.D-c.U) \ (c.L*Px) + h;
%         r = Pf - PAP * Px;
        
        itr = itr + 1;
        if mod(itr, 100) == 0
            fprintf('Iteration %d / %.3e.\n', itr, norm(r)/norm(r0));
        end
        
        if itr >= itrLimit || (norm(r) < norm(r0) * 1e-6)
            break
        end
    end
    
    x(c.dul) = Px(c.dlu);
    x = x';
%     x = c.P * Px;
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

function M = createTriDiagSparse(m, q, n)
    i = [1:n 1:n-1 2:n];
    j = [1:n 2:n 1:n-1];
    v = [ones(1, n)*m ones(1, 2*(n-1))*q];
    M = sparse(i, j, v);
end